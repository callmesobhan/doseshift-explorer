"""Minimal synthetic DICOM-RT writer for tests and demos. Not conformant enough for clinical use."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pydicom
from contourpy import contour_generator
from pydicom.dataset import Dataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

from ..synthetic import PhantomCase

CT_SOP = "1.2.840.10008.5.1.4.1.1.2"
RTDOSE_SOP = "1.2.840.10008.5.1.4.1.1.481.2"
RTSTRUCT_SOP = "1.2.840.10008.5.1.4.1.1.481.3"
RTPLAN_SOP = "1.2.840.10008.5.1.4.1.1.481.5"


def _base(sop_class: str, modality: str, uids: dict[str, str]) -> Dataset:
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = sop_class
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian
    meta.ImplementationClassUID = pydicom.uid.PYDICOM_IMPLEMENTATION_UID
    ds = Dataset()
    ds.file_meta = meta
    ds.SOPClassUID = sop_class
    ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
    ds.Modality = modality
    ds.PatientName = "DoseShift^Phantom"
    ds.PatientID = "PHANTOM"
    ds.StudyInstanceUID = uids["study"]
    ds.FrameOfReferenceUID = uids["for"]
    return ds


def _geometry(ds: Dataset, grid, k_index_origin: int = 0) -> None:
    a = grid.affine
    spacing = grid.spacing_xyz
    ds.ImageOrientationPatient = [*(a[:3, 0] / spacing[0]), *(a[:3, 1] / spacing[1])]
    ds.PixelSpacing = [float(spacing[1]), float(spacing[0])]  # [row, column]
    ds.Rows, ds.Columns = grid.shape[1], grid.shape[2]


def _set_pixels(ds: Dataset, array: np.ndarray, vr: str = "OW") -> None:
    ds.PixelData = array.tobytes()
    ds["PixelData"].VR = vr


def write_ct(case: PhantomCase, outdir: Path, uids: dict[str, str]) -> Path:
    out = outdir / "CT"
    out.mkdir(parents=True, exist_ok=True)
    g = case.ct_grid
    series = uids["ct_series"] = generate_uid()
    for k in range(g.shape[0]):
        ds = _base(CT_SOP, "CT", uids)
        ds.SeriesInstanceUID = series
        ds.InstanceNumber = k + 1
        _geometry(ds, g)
        ds.ImagePositionPatient = [float(v) for v in g.index_to_world([[0, 0, k]])[0]]
        ds.SliceThickness = float(g.spacing_xyz[2])
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.BitsAllocated, ds.BitsStored, ds.HighBit, ds.PixelRepresentation = 16, 16, 15, 0
        ds.RescaleIntercept, ds.RescaleSlope = -1024.0, 1.0
        raw = np.clip(np.rint(case.ct_hu[k] + 1024.0), 0, 65535).astype("<u2")
        _set_pixels(ds, raw)
        ds.save_as(out / f"CT_{k:03d}.dcm", enforce_file_format=True)
    return out


def write_rtdose(
    case: PhantomCase, outdir: Path, uids: dict[str, str], reverse_z: bool = False
) -> Path:
    """``reverse_z`` stores frames from +z to -z (negative GridFrameOffsetVector steps)."""
    g, data = case.dose.grid, case.dose.data
    nz = g.shape[0]
    ds = _base(RTDOSE_SOP, "RTDOSE", uids)
    ds.SeriesInstanceUID = generate_uid()
    ds.InstanceNumber = 1
    _geometry(ds, g)
    z0 = g.spacing_xyz[2]
    if reverse_z:
        data = data[::-1]
        first = g.index_to_world([[0, 0, nz - 1]])[0]
        offsets = [-k * z0 for k in range(nz)]
    else:
        first = g.index_to_world([[0, 0, 0]])[0]
        offsets = [k * z0 for k in range(nz)]
    ds.ImagePositionPatient = [float(v) for v in first]
    ds.GridFrameOffsetVector = [float(o) for o in offsets]
    ds.FrameIncrementPointer = 0x3004000C
    ds.NumberOfFrames = nz
    if "plan_sop" in uids:
        ref = Dataset()
        ref.ReferencedSOPClassUID, ref.ReferencedSOPInstanceUID = RTPLAN_SOP, uids["plan_sop"]
        ds.ReferencedRTPlanSequence = Sequence([ref])
    ds.DoseUnits, ds.DoseType, ds.DoseSummationType = "GY", "PHYSICAL", "PLAN"
    scaling = float(data.max()) / (2**32 - 1)
    ds.DoseGridScaling = scaling
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated, ds.BitsStored, ds.HighBit, ds.PixelRepresentation = 32, 32, 31, 0
    _set_pixels(ds, np.rint(data / scaling).astype("<u4"))
    path = outdir / ("RTDOSE_reversed.dcm" if reverse_z else "RTDOSE.dcm")
    ds.save_as(path, enforce_file_format=True)
    return path


def _mask_contours(mask2d: np.ndarray) -> list[np.ndarray]:
    """Closed polygons (n,2) in (i, j) index space tracing the 0.5 level of a binary slice."""
    padded = np.pad(mask2d.astype(float), 1)
    polys = []
    for line in contour_generator(z=padded).lines(0.5):
        pts = np.asarray(line) - 1.0  # undo the padding
        polys.append(pts[:-1] if np.allclose(pts[0], pts[-1]) else pts)
    return polys


def write_rtstruct(case: PhantomCase, outdir: Path, uids: dict[str, str]) -> Path:
    ds = _base(RTSTRUCT_SOP, "RTSTRUCT", uids)
    ds.SeriesInstanceUID = generate_uid()
    ds.StructureSetLabel = "PHANTOM"
    roi_seq, contour_seq, obs_seq = Sequence(), Sequence(), Sequence()
    for n, (name, s) in enumerate(case.masks.items(), start=1):
        roi = Dataset()
        roi.ROINumber, roi.ROIName = n, name
        roi.ReferencedFrameOfReferenceUID = uids["for"]
        roi.ROIGenerationAlgorithm = "MANUAL"
        roi_seq.append(roi)
        rc = Dataset()
        rc.ReferencedROINumber = n
        rc.ROIDisplayColor = [255, 0, 0]
        rc.ContourSequence = Sequence()
        for k in range(s.mask.shape[0]):
            if not s.mask[k].any():
                continue
            for poly in _mask_contours(s.mask[k]):
                ijk = np.column_stack([poly, np.full(len(poly), k)])
                c = Dataset()
                c.ContourGeometricType = "CLOSED_PLANAR"
                c.NumberOfContourPoints = len(poly)
                c.ContourData = [float(v) for v in s.grid.index_to_world(ijk).ravel()]
                rc.ContourSequence.append(c)
        contour_seq.append(rc)
        ob = Dataset()
        ob.ObservationNumber, ob.ReferencedROINumber = n, n
        ob.RTROIInterpretedType = {"ptv": "PTV", "ctv": "CTV"}.get(s.role, "ORGAN")
        obs_seq.append(ob)
    series_ref = Dataset()
    series_ref.SeriesInstanceUID = uids["ct_series"]
    study_ref = Dataset()
    study_ref.ReferencedSOPClassUID = "1.2.840.10008.3.1.2.3.1"
    study_ref.ReferencedSOPInstanceUID = uids["study"]
    study_ref.RTReferencedSeriesSequence = Sequence([series_ref])
    frame_ref = Dataset()
    frame_ref.FrameOfReferenceUID = uids["for"]
    frame_ref.RTReferencedStudySequence = Sequence([study_ref])
    ds.ReferencedFrameOfReferenceSequence = Sequence([frame_ref])
    uids["struct_sop"] = str(ds.SOPInstanceUID)
    ds.StructureSetROISequence = roi_seq
    ds.ROIContourSequence = contour_seq
    ds.RTROIObservationsSequence = obs_seq
    path = outdir / "RTSTRUCT.dcm"
    ds.save_as(path, enforce_file_format=True)
    return path


def write_rtplan(
    case: PhantomCase, outdir: Path, uids: dict[str, str], fractions: int = 20
) -> Path:
    """Minimal RTPLAN: prescription, fractions, and the structure set it was planned on."""
    ds = _base(RTPLAN_SOP, "RTPLAN", uids)
    ds.SeriesInstanceUID = generate_uid()
    ds.RTPlanLabel = "PHANTOM"
    ref = Dataset()
    ref.DoseReferenceNumber, ref.DoseReferenceType = 1, "TARGET"
    ref.TargetPrescriptionDose = float(case.prescription_gy)
    ds.DoseReferenceSequence = Sequence([ref])
    fg = Dataset()
    fg.FractionGroupNumber, fg.NumberOfFractionsPlanned = 1, fractions
    ds.FractionGroupSequence = Sequence([fg])
    sref = Dataset()
    sref.ReferencedSOPClassUID, sref.ReferencedSOPInstanceUID = RTSTRUCT_SOP, uids["struct_sop"]
    ds.ReferencedStructureSetSequence = Sequence([sref])
    uids["plan_sop"] = str(ds.SOPInstanceUID)
    path = outdir / "RTPLAN.dcm"
    ds.save_as(path, enforce_file_format=True)
    return path


def write_case(
    case: PhantomCase, outdir: str | Path, reverse_dose_z: bool = False
) -> dict[str, Path]:
    """Write CT/, RTDOSE.dcm and RTSTRUCT.dcm for a synthetic phantom."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    uids = {"study": generate_uid(), "for": generate_uid()}
    ct_dir = write_ct(case, outdir, uids)
    rtstruct = write_rtstruct(case, outdir, uids)
    rtplan = write_rtplan(case, outdir, uids)
    dose = write_rtdose(case, outdir, uids, reverse_z=reverse_dose_z)
    return {"ct_dir": ct_dir, "dose": dose, "rtstruct": rtstruct, "rtplan": rtplan}
