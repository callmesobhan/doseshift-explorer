"""Read CT, RTDOSE and RTSTRUCT into the project's Grid / mask / PlanState types.

Geometry follows DICOM PS3.3 C.7.6.2: ImagePositionPatient is the *centre* of the first voxel,
ImageOrientationPatient gives the row (x) and column (y) direction cosines, PixelSpacing is
[row spacing, column spacing]. RTDOSE slice positions come from GridFrameOffsetVector (offsets from
ImagePositionPatient), which may run against the slice normal; the sign is preserved in the affine.

Malformed input raises ``ValueError`` (or ``FileNotFoundError``) naming the file and the problem;
suspicious but usable input (dropped contours, mismatched references) warns.
"""

from __future__ import annotations

import warnings
from dataclasses import replace
from pathlib import Path

import numpy as np
import pydicom
from matplotlib.path import Path as MplPath

from ..geometry import Grid
from ..state import DoseGrid, MaskStructure, PlanState

_PLANAR = ("CLOSED_PLANAR", "INTERPOLATED_PLANAR")
IMAGE_MODALITIES = ("CT", "MR")


def _require(ds, tags: tuple[str, ...], path) -> None:
    missing = [t for t in tags if t not in ds]
    if missing:
        raise ValueError(f"{path}: missing required DICOM attribute(s) {missing}")


def _pixels(ds, path) -> np.ndarray:
    try:
        return ds.pixel_array
    except Exception as exc:  # pydicom raises many types depending on the missing decoder
        ts = getattr(getattr(ds, "file_meta", None), "TransferSyntaxUID", "unknown")
        raise ValueError(
            f"{path}: cannot decode pixel data (transfer syntax {ts}): {exc}. Compressed images "
            "need an extra decoder, e.g. `pip install pylibjpeg pylibjpeg-libjpeg`"
        ) from exc


def _iop_axes(iop, path="") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    iop = np.asarray(iop, float)
    if iop.shape != (6,) or not np.isfinite(iop).all():
        raise ValueError(f"{path}: invalid ImageOrientationPatient {iop.tolist()}")
    row_dir, col_dir = iop[:3], iop[3:]
    normal = np.cross(row_dir, col_dir)
    if not np.isclose(np.linalg.norm(normal), 1.0, atol=1e-2):
        raise ValueError(f"{path}: ImageOrientationPatient is not two orthonormal vectors")
    return row_dir, col_dir, normal


def _affine(origin, row_dir, col_dir, k_vec, pixel_spacing, path="") -> np.ndarray:
    ps = np.asarray(pixel_spacing, float)
    if ps.shape != (2,) or not (np.isfinite(ps).all() and (ps > 0).all()):
        raise ValueError(f"{path}: invalid PixelSpacing {ps.tolist()}")
    a = np.eye(4)
    a[:3, 0] = row_dir * ps[1]  # i runs along the row direction, step = column spacing
    a[:3, 1] = col_dir * ps[0]  # j runs along the column direction, step = row spacing
    a[:3, 2] = k_vec
    a[:3, 3] = np.asarray(origin, float)
    return a


def load_image_series(
    directory: str | Path, modalities: tuple[str, ...] = IMAGE_MODALITIES
) -> tuple[Grid, np.ndarray, str]:
    """CT or MR series -> (grid, values [k, j, i], FrameOfReferenceUID).

    Values are HU for CT (rescale applied) and raw scaled intensities for MR. The directory must
    hold exactly one series with uniform slice spacing; non-DICOM files are ignored.
    """
    directory = Path(directory)
    if not directory.is_dir():
        raise FileNotFoundError(f"image directory not found: {directory}")
    slices = []
    for f in sorted(directory.iterdir()):
        if f.is_file():
            try:
                ds = pydicom.dcmread(f)
            except pydicom.errors.InvalidDicomError:
                continue
            if ds.get("Modality") in modalities and "ImagePositionPatient" in ds:
                ds.filename = str(f)
                slices.append(ds)
    if not slices:
        raise FileNotFoundError(f"no {'/'.join(modalities)} slices found in {directory}")
    for ds in slices:
        _require(
            ds,
            ("ImageOrientationPatient", "PixelSpacing", "SeriesInstanceUID", "Rows", "Columns"),
            ds.filename,
        )
    if len({str(d.SeriesInstanceUID) for d in slices}) > 1:
        raise ValueError(f"{directory} holds more than one series; give one series per directory")
    if len({(int(d.Rows), int(d.Columns)) for d in slices}) > 1:
        raise ValueError(f"{directory}: slices have different Rows/Columns")
    iops = np.array([np.asarray(d.ImageOrientationPatient, float) for d in slices])
    if not np.allclose(iops, iops[0], atol=1e-3):
        raise ValueError(f"{directory}: slices have different ImageOrientationPatient")
    row_dir, col_dir, normal = _iop_axes(iops[0], directory)
    slices.sort(key=lambda d: float(np.dot(np.asarray(d.ImagePositionPatient, float), normal)))
    pos = np.array([np.asarray(d.ImagePositionPatient, float) for d in slices])
    if not np.isfinite(pos).all():
        raise ValueError(f"{directory}: non-finite ImagePositionPatient")
    if len(slices) > 1:
        steps = np.diff(pos @ normal)
        if not np.allclose(steps, steps[0], atol=1e-2) or abs(steps[0]) < 1e-6:
            raise ValueError("non-uniform (or zero) slice spacing is not supported")
        k_vec = (pos[-1] - pos[0]) / (len(slices) - 1)
    else:
        k_vec = normal * float(slices[0].get("SliceThickness", 1.0))
    affine = _affine(pos[0], row_dir, col_dir, k_vec, slices[0].PixelSpacing, directory)
    values = np.stack(
        [
            _pixels(d, d.filename).astype(np.float32) * float(d.get("RescaleSlope", 1))
            + float(d.get("RescaleIntercept", 0))
            for d in slices
        ]
    )
    return Grid(values.shape, affine), values, str(slices[0].get("FrameOfReferenceUID", ""))


def load_ct(directory: str | Path) -> tuple[Grid, np.ndarray, str]:
    """CT series -> (grid, HU, FrameOfReferenceUID); see ``load_image_series``."""
    return load_image_series(directory, ("CT",))


def series_uid(directory: str | Path) -> str:
    """SeriesInstanceUID of the (single) series stored in ``directory``."""
    for f in sorted(Path(directory).iterdir()):
        if f.is_file():
            try:
                return str(pydicom.dcmread(f, stop_before_pixels=True).SeriesInstanceUID)
            except (pydicom.errors.InvalidDicomError, AttributeError):
                continue
    raise FileNotFoundError(f"no DICOM files in {directory}")


def rtstruct_referenced_series(ds) -> list[str]:
    """Image SeriesInstanceUIDs an RTSTRUCT says its contours were drawn on."""
    out = []
    for fr in ds.get("ReferencedFrameOfReferenceSequence", []):
        for st in fr.get("RTReferencedStudySequence", []):
            out += [str(se.SeriesInstanceUID) for se in st.get("RTReferencedSeriesSequence", [])]
    return out


def check_references(image_dir: str | Path, rtstruct_path: str | Path) -> str | None:
    """Warning text if the RTSTRUCT does not reference the image series in ``image_dir``."""
    refs = rtstruct_referenced_series(pydicom.dcmread(rtstruct_path, stop_before_pixels=True))
    uid = series_uid(image_dir)
    if refs and uid not in refs:
        return (
            f"RTSTRUCT references image series ...{[r[-8:] for r in refs]} but the loaded image "
            f"series is ...{uid[-8:]}; contours may not belong to this image"
        )
    return None


def load_rtdose(path: str | Path) -> tuple[DoseGrid, str]:
    """RTDOSE -> (DoseGrid in Gy, FrameOfReferenceUID)."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"RTDOSE file not found: {path}")
    ds = pydicom.dcmread(path)
    _require(
        ds,
        ("DoseGridScaling", "ImagePositionPatient", "ImageOrientationPatient", "PixelSpacing"),
        path,
    )
    if str(ds.get("DoseUnits", "GY")).upper() != "GY":
        raise ValueError(f"{path}: DoseUnits={ds.DoseUnits!r}; only absolute Gy is supported")
    scaling = float(ds.DoseGridScaling)
    if not (np.isfinite(scaling) and scaling > 0):
        raise ValueError(f"{path}: DoseGridScaling must be finite and > 0, got {scaling}")
    data = np.asarray(_pixels(ds, path), dtype=float)
    if data.ndim == 2:
        data = data[None]
    data = data * scaling
    if not np.isfinite(data).all() or data.min() < 0:
        raise ValueError(f"{path}: dose contains non-finite or negative values")
    row_dir, col_dir, normal = _iop_axes(ds.ImageOrientationPatient, path)
    offsets = np.asarray(ds.get("GridFrameOffsetVector", [0.0]), dtype=float)
    if len(offsets) != data.shape[0]:
        raise ValueError(f"{path}: GridFrameOffsetVector length does not match number of frames")
    if len(offsets) > 1:
        steps = np.diff(offsets)
        if not np.allclose(steps, steps[0], atol=1e-3) or abs(steps[0]) < 1e-6:
            raise ValueError(f"{path}: non-uniform GridFrameOffsetVector is not supported")
        dz = steps[0]  # signed along the slice normal
    else:
        dz = float(ds.get("SliceThickness", 1.0))
    # offsets are measured from ImagePositionPatient; the first frame is usually at offset 0
    origin = np.asarray(ds.ImagePositionPatient, float) + normal * offsets[0]
    affine = _affine(origin, row_dir, col_dir, normal * dz, ds.PixelSpacing, path)
    return DoseGrid(Grid(data.shape, affine), data), str(ds.get("FrameOfReferenceUID", ""))


def rasterize_contours(contours_world: list[np.ndarray], grid: Grid, label: str = "") -> np.ndarray:
    """Planar contours (list of (n,3) world polygons) -> bool mask. A voxel is inside if its
    centre is inside the polygon; overlapping contours on a slice combine by XOR (holes)."""
    nz, ny, nx = grid.shape
    tag = f"{label}: " if label else ""
    mask = np.zeros(grid.shape, bool)
    skipped = 0
    for pts in contours_world:
        if not np.isfinite(pts).all():
            raise ValueError(f"{tag}contour has non-finite (NaN/inf) coordinates")
        ijk = grid.world_to_index(pts)
        k_float = ijk[:, 2].mean()
        k = int(round(k_float))
        if not (0 <= k < nz) or abs(k_float - k) > 0.25:
            skipped += 1
            continue
        i0, j0 = np.floor(ijk[:, :2].min(axis=0)).astype(int)
        i1, j1 = np.ceil(ijk[:, :2].max(axis=0)).astype(int)
        i0, j0, i1, j1 = max(i0, 0), max(j0, 0), min(i1, nx - 1), min(j1, ny - 1)
        if i1 < i0 or j1 < j0:
            continue
        jj, ii = np.meshgrid(np.arange(j0, j1 + 1), np.arange(i0, i1 + 1), indexing="ij")
        inside = MplPath(ijk[:, :2]).contains_points(np.stack([ii.ravel(), jj.ravel()], 1))
        mask[k, jj.ravel(), ii.ravel()] ^= inside
    if skipped:
        warnings.warn(
            f"{tag}{skipped} of {len(contours_world)} contour(s) did not lie on a slice of the "
            "image grid and were skipped (wrong image series for this RTSTRUCT?)",
            stacklevel=2,
        )
    return mask


def load_rtstruct_masks(
    path: str | Path, grid: Grid, names: list[str] | None = None
) -> tuple[dict[str, np.ndarray], str]:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"RTSTRUCT file not found: {path}")
    ds = pydicom.dcmread(path)
    _require(ds, ("StructureSetROISequence", "ROIContourSequence"), path)
    roi_names: dict[int, str] = {}
    for r in ds.StructureSetROISequence:
        roi_names[int(r.ROINumber)] = str(r.ROIName)
    wanted = None if names is None else set(names)
    seen: dict[str, int] = {}
    masks: dict[str, np.ndarray] = {}
    for rc in ds.ROIContourSequence:
        number = int(rc.ReferencedROINumber)
        if number not in roi_names:
            raise ValueError(f"{path}: contour references unknown ROI number {number}")
        name = roi_names[number]
        if wanted is not None and name not in wanted:
            continue
        seen[name] = seen.get(name, 0) + 1
        if seen[name] > 1:
            raise ValueError(
                f"{path}: more than one ROI is named {name!r}; rename one in the planning system "
                "(masks cannot be told apart by name)"
            )
        polys, dropped = [], 0
        for c in getattr(rc, "ContourSequence", []):
            if str(c.ContourGeometricType) not in _PLANAR:
                dropped += 1
                continue
            data = np.asarray(c.ContourData, float)
            if data.size % 3 or data.size < 9:
                raise ValueError(f"{name!r}: malformed ContourData (length {data.size})")
            polys.append(data.reshape(-1, 3))
        if dropped:
            warnings.warn(
                f"{name}: {dropped} non-planar contour(s) (e.g. POINT/OPEN_PLANAR) ignored",
                stacklevel=2,
            )
        masks[name] = rasterize_contours(polys, grid, name)
    for_uid = str(ds.get("FrameOfReferenceUID", "")) or _struct_for_uid(ds)
    return masks, for_uid


def _struct_for_uid(ds) -> str:
    seq = ds.get("ReferencedFrameOfReferenceSequence")
    return str(seq[0].FrameOfReferenceUID) if seq else ""


def load_plan_state_with_image(
    image_dir: str | Path,
    dose_path: str | Path,
    rtstruct_path: str | Path,
    roles: dict[str, str],
    prescription_gy: float,
    rx_label: str = "Rx",
) -> tuple[PlanState, tuple[Grid, np.ndarray]]:
    """Like ``load_plan_state`` but also returns the image ``(grid, values)``, loaded once."""
    ct_grid, values, ct_for = load_image_series(image_dir)
    dose, dose_for = load_rtdose(dose_path)
    masks, rs_for = load_rtstruct_masks(rtstruct_path, ct_grid, list(roles))
    missing = set(roles) - set(masks)
    if missing:
        raise KeyError(f"ROIs not found in RTSTRUCT: {sorted(missing)}")
    mismatch = check_references(image_dir, rtstruct_path)
    if mismatch:
        warnings.warn(mismatch, stacklevel=2)
    for label, uid in (("RTDOSE", dose_for), ("RTSTRUCT", rs_for)):
        if uid and ct_for and uid != ct_for:
            warnings.warn(f"{label} FrameOfReferenceUID differs from the CT's", stacklevel=2)
    structures = [MaskStructure(n, roles[n], masks[n], ct_grid) for n in roles]
    state = PlanState.from_masks(dose, structures, prescription_gy)
    return replace(state, rx_label=rx_label), (ct_grid, values)


def load_plan_state(
    ct_dir: str | Path,
    dose_path: str | Path,
    rtstruct_path: str | Path,
    roles: dict[str, str],
    prescription_gy: float,
    rx_label: str = "Rx",
) -> PlanState:
    """Build a PlanState; ``roles`` maps ROI name -> 'ptv' | 'ctv' | 'oar' (others are ignored)."""
    state, _ = load_plan_state_with_image(
        ct_dir, dose_path, rtstruct_path, roles, prescription_gy, rx_label
    )
    return state


def alignment_report(
    ct_dir: str | Path, dose_path: str | Path, rtstruct_path: str | Path | None = None
) -> list[str]:
    """Human-readable checks of the dose-to-CT geometric relationship."""
    ct_grid, _, ct_for = load_image_series(ct_dir)
    dose, dose_for = load_rtdose(dose_path)
    g = dose.grid
    lines = [
        f"CT   grid {ct_grid.shape} spacing {np.round(ct_grid.spacing_xyz, 3)} mm",
        f"DOSE grid {g.shape} spacing {np.round(g.spacing_xyz, 3)} mm",
        f"FrameOfReference: CT == DOSE : {ct_for == dose_for}",
    ]
    if rtstruct_path is not None:
        _, rs_for = load_rtstruct_masks(rtstruct_path, ct_grid, names=[])
        lines.append(f"FrameOfReference: CT == RTSTRUCT : {ct_for == rs_for}")
        lines.append(
            "RTSTRUCT references the loaded image series : "
            + ("yes" if check_references(ct_dir, rtstruct_path) is None else "NO (see warning)")
        )
    phase = ct_grid.world_to_index(g.origin_xyz[None])[0]
    lines.append(
        "dose voxel (0,0,0) sits at CT index "
        f"{np.round(phase, 3)} (fractional part = sub-voxel offset between the two grids)"
    )
    lo, hi = g.extent_xyz()
    lines.append(
        f"dose extent (voxel centres) x[{lo[0]:.1f},{hi[0]:.1f}] "
        f"y[{lo[1]:.1f},{hi[1]:.1f}] z[{lo[2]:.1f},{hi[2]:.1f}] mm"
    )
    return lines
