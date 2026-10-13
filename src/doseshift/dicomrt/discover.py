"""Assemble analysable cases from a folder of downloaded DICOM series.

Downloads such as TCIA collections give one directory per series with nothing in the file names
saying which RTSTRUCT belongs to which image or which RTDOSE to which plan. This module follows
the references stored in the DICOM headers instead of guessing from names:

    RTDOSE --ReferencedRTPlan--> RTPLAN --ReferencedStructureSet--> RTSTRUCT
                                                                     --ReferencedSeries--> image

When no RTPLAN is present (some collections omit it), an RTDOSE is paired by FrameOfReference only
if that is unambiguous; otherwise the case is reported without a dose and the reason is listed.
Only headers are read (no pixel data).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pydicom

from .reader import rtstruct_referenced_series

IMAGE_MODALITIES = ("CT", "MR")


@dataclass
class CaseFiles:
    patient: str
    label: str
    image_dir: Path
    rtstruct: Path
    rtdose: Path | None
    rtplan: Path | None
    prescription_gy: float | None
    n_fractions: int | None
    roi_names: list[str]
    notes: list[str] = field(default_factory=list)

    @property
    def complete(self) -> bool:
        return self.rtdose is not None

    def summary(self) -> str:
        rx = f"{self.prescription_gy:g} Gy" if self.prescription_gy is not None else "unknown"
        dose = self.rtdose.parent.name if self.rtdose else "NONE"
        lines = [
            f"[{self.patient}] {self.label}",
            f"    image   {self.image_dir.name}",
            f"    rtstruct {self.rtstruct.parent.name}  ROIs: {', '.join(self.roi_names)}",
            f"    rtdose  {dose}   rtplan {self.rtplan.parent.name if self.rtplan else '-'}"
            f"   prescription {rx}" + (f" in {self.n_fractions} fx" if self.n_fractions else ""),
        ]
        lines += [f"    note: {n}" for n in self.notes]
        return "\n".join(lines)


def _header(path: Path):
    return pydicom.dcmread(path, stop_before_pixels=True)


def _dicom_files(root: Path) -> list[Path]:
    """All DICOM files below ``root`` whatever their extension (TPS exports often have none)."""
    out = []
    for f in sorted(root.rglob("*")):
        if f.is_file() and not f.is_symlink() and pydicom.misc.is_dicom(f):
            out.append(f)
    return out


def _ref_uid(ds, seq_name: str) -> str | None:
    seq = ds.get(seq_name)
    return str(seq[0].ReferencedSOPInstanceUID) if seq else None


def _prescription(plan) -> tuple[float | None, int | None]:
    rx = None
    for ref in plan.get("DoseReferenceSequence", []):
        value = ref.get("TargetPrescriptionDose")
        if value is not None and str(ref.get("DoseReferenceType", "TARGET")).upper() == "TARGET":
            rx = max(rx or 0.0, float(value))
    fx = None
    for fg in plan.get("FractionGroupSequence", []):
        n = fg.get("NumberOfFractionsPlanned")
        fx = int(n) if n is not None and fx is None else fx
    return rx, fx


def discover_cases(patient_dir: str | Path) -> list[CaseFiles]:
    root = Path(patient_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"directory not found: {root}")
    # one entry per image series (the directory) and one per RT object (the file), so both
    # one-series-per-folder layouts (TCIA) and flat folders of RT files work
    entries: list[tuple[Path, object]] = []
    by_dir: dict[Path, list[Path]] = {}
    for f in _dicom_files(root):
        by_dir.setdefault(f.parent, []).append(f)
    for d, files in by_dir.items():
        first = _header(files[0])
        if first.Modality in IMAGE_MODALITIES:
            entries.append((d, first))
        else:
            entries += [(files[0], first)] + [(f, _header(f)) for f in files[1:]]

    images = {str(ds.SeriesInstanceUID): d for d, ds in entries if ds.Modality in IMAGE_MODALITIES}
    plans = {str(ds.SOPInstanceUID): (f, ds) for f, ds in entries if ds.Modality == "RTPLAN"}
    doses = [(f, ds) for f, ds in entries if ds.Modality == "RTDOSE"]
    patient = root.name

    cases = []
    for f, rs in ((f, ds) for f, ds in entries if ds.Modality == "RTSTRUCT"):
        notes: list[str] = []
        refs = rtstruct_referenced_series(rs)
        image_dir = next((images[r] for r in refs if r in images), None)
        if image_dir is None:
            continue  # contours for an image series that was not downloaded
        image_ds = dict(entries)[image_dir]
        rois = [str(r.ROIName) for r in rs.get("StructureSetROISequence", [])]
        struct_uid = str(rs.SOPInstanceUID)

        plan_f = plan_ds = None
        for pf, pds in plans.values():
            if _ref_uid(pds, "ReferencedStructureSetSequence") == struct_uid:
                plan_f, plan_ds = pf, pds
                break
        dose_f = None
        if plan_ds is not None:
            for df, dds in doses:
                if _ref_uid(dds, "ReferencedRTPlanSequence") == str(plan_ds.SOPInstanceUID):
                    dose_f = df
                    break
            if dose_f is None:
                notes.append("RTPLAN found but no RTDOSE references it")
        else:
            same_for = [
                df
                for df, dds in doses
                if str(dds.get("FrameOfReferenceUID"))
                == str(rs.get("FrameOfReferenceUID", image_ds.get("FrameOfReferenceUID")))
            ]
            if len(same_for) == 1:
                dose_f = same_for[0]
                notes.append("no RTPLAN; RTDOSE paired by FrameOfReferenceUID (unique match)")
            elif same_for:
                notes.append(
                    f"{len(same_for)} RTDOSE share the FrameOfReference and no RTPLAN "
                    "disambiguates them; not paired"
                )
            else:
                notes.append("no matching RTDOSE")
        rx, fx = _prescription(plan_ds) if plan_ds is not None else (None, None)
        if rx is None:
            notes.append("prescription not in the data: supply prescription_gy yourself")
        desc = str(image_ds.get("SeriesDescription", "") or image_ds.Modality)
        cases.append(
            CaseFiles(
                patient,
                f"{rs.get('StructureSetLabel', '')} on {image_ds.Modality} '{desc}'",
                image_dir,
                f,
                dose_f,
                plan_f,
                rx,
                fx,
                rois,
                notes,
            )
        )
    return cases


_ABSOLUTE_METRICS = {
    "ptv": ["D98%", "D95%", "D50%", "D2%", "Dmean", "HI"],
    "ctv": ["D98%", "D95%", "Dmean"],
    "oar": ["Dmean", "D2%", "D5cc", "D0.5cc"],
}


def case_config(
    case: CaseFiles,
    roles: dict[str, str],
    experiment: str = "grid_audit",
    prescription_gy: float | None = None,
) -> dict:
    """A ready-to-run experiment config for ``case``; ``roles`` maps ROI name -> ptv|ctv|oar."""
    if not case.complete:
        raise ValueError("case has no RTDOSE: " + "; ".join(case.notes))
    rx = prescription_gy if prescription_gy is not None else case.prescription_gy
    if rx is None:
        raise ValueError("prescription unknown: pass prescription_gy")
    unknown = set(roles) - set(case.roi_names)
    if unknown:
        raise ValueError(f"ROIs not in this RTSTRUCT: {sorted(unknown)}; has {case.roi_names}")
    ptv = next((n for n, r in roles.items() if r == "ptv"), None)
    oar = next((n for n, r in roles.items() if r == "oar"), None)
    keys = [f"{ptv}:D95%"] if ptv else []
    keys += [f"{oar}:D2%"] if oar else []
    cfg: dict = {
        "experiment": experiment,
        "case": {
            "type": "dicom",
            "image_dir": str(case.image_dir),
            "dose": str(case.rtdose),
            "rtstruct": str(case.rtstruct),
            "prescription_gy": float(rx),
            "roles": dict(roles),
        },
    }
    if case.prescription_gy is None:
        cfg["case"]["prescription_is_reference"] = True
        # a user-supplied reference dose is not a prescription: drop metrics that depend on it
        cfg["metrics"] = {n: list(_ABSOLUTE_METRICS[r]) for n, r in roles.items()}
    if experiment == "grid_audit":
        cfg["audit"] = {"keys": keys}
    elif experiment == "sweep":
        cfg["sweep"] = {"keys": keys}
    return cfg
