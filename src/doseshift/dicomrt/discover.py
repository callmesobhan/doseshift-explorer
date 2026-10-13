"""Assemble analysable cases from a folder of downloaded DICOM series.

Downloads such as TCIA collections give one directory per series with nothing in the file names
saying which RTSTRUCT belongs to which image or which RTDOSE to which plan. This module follows
the references stored in the DICOM headers instead of guessing from names:

    RTDOSE --ReferencedRTPlan--> RTPLAN --ReferencedStructureSet--> RTSTRUCT
                                                                     --ReferencedSeries--> image

Ambiguity is never resolved silently:

* one case per (RTSTRUCT, image series, RTPLAN) that is actually present;
* among several RTDOSE of one plan, plan-level physical dose wins over fraction/beam/effective
  dose; if candidates are still tied the case gets no dose and says why; ignored doses are listed;
* a prescription is reported only when the plan states exactly one target dose in one phase;
* RTSTRUCTs that cannot be placed (no image reference, image not in the folder) are returned in a
  report instead of vanishing (``discover_with_report``; ``discover_cases`` warns about them).

When no RTPLAN references the structure set, an RTDOSE is paired by FrameOfReference only if that
is unambiguous. Images are grouped by SeriesInstanceUID, not by folder. Only headers are read.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import pydicom

from .reader import rtstruct_referenced_series

IMAGE_MODALITIES = ("CT", "MR")
_PLAN_LEVEL = ("PLAN", "MULTI_PLAN")


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
    """All DICOM files below ``root`` whatever their extension (TPS exports often have none).
    Symlinked files and directories are followed (each real directory is visited once)."""
    out, seen = [], set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen:
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames.sort()
        for name in sorted(filenames):
            f = Path(dirpath) / name
            if f.is_file() and pydicom.misc.is_dicom(f):
                out.append(f)
    return out


def _ref_uids(ds, seq_name: str) -> list[str]:
    """Every ReferencedSOPInstanceUID in a sequence (an object may reference several)."""
    return [str(item.ReferencedSOPInstanceUID) for item in ds.get(seq_name, [])]


def _dedupe(items):
    """Drop files that repeat an already seen SOPInstanceUID; returns (unique, {uid: [names]})."""
    unique, dups = {}, {}
    for f, ds in items:
        uid = str(ds.SOPInstanceUID)
        if uid in unique:
            dups.setdefault(uid, []).append(f.name)
        else:
            unique[uid] = (f, ds)
    return list(unique.values()), dups


def _short(uid: str) -> str:
    return "..." + uid[-8:]


def _prescription(plan, notes: list[str]) -> tuple[float | None, int | None]:
    """(prescription Gy, fractions): only when the plan states one target dose in one phase."""
    targets = sorted(
        {
            float(ref.TargetPrescriptionDose)
            for ref in plan.get("DoseReferenceSequence", [])
            if ref.get("TargetPrescriptionDose") is not None
            and str(ref.get("DoseReferenceType", "TARGET")).upper() == "TARGET"
        }
    )
    groups = list(plan.get("FractionGroupSequence", []))
    fractions = [
        int(g.NumberOfFractionsPlanned)
        for g in groups
        if g.get("NumberOfFractionsPlanned") is not None
    ]
    n_fx = sum(fractions) if fractions else None
    if len(groups) > 1:
        notes.append(
            f"multi-phase plan ({len(groups)} fraction groups, fractions summed): the prescription "
            "differs per phase, supply prescription_gy yourself"
        )
        return None, n_fx
    if len(targets) > 1:
        notes.append(
            f"plan states {len(targets)} different target doses {targets} Gy: "
            "supply prescription_gy"
        )
        return None, n_fx
    if targets and targets[0] <= 0:
        notes.append(f"plan prescription is {targets[0]} Gy (invalid)")
        return None, n_fx
    return (targets[0] if targets else None), n_fx


def _rank(ds) -> int:
    """Higher is the better dose to analyse against a plan prescription."""
    summation = str(ds.get("DoseSummationType", "")).upper()
    score = 2 if summation in _PLAN_LEVEL else 1 if summation.startswith("FRACTION") else 0
    return 2 * score + (1 if str(ds.get("DoseType", "PHYSICAL")).upper() == "PHYSICAL" else 0)


def _choose_dose(candidates: list[tuple[Path, object]], notes: list[str]) -> Path | None:
    if not candidates:
        return None
    best = max(_rank(ds) for _, ds in candidates)
    top = [(f, ds) for f, ds in candidates if _rank(ds) == best]
    rest = [(f, ds) for f, ds in candidates if _rank(ds) != best]
    if rest:
        kinds = sorted({str(ds.get("DoseSummationType", "?")) for _, ds in rest})
        notes.append(f"ignored {len(rest)} lower-priority RTDOSE ({', '.join(kinds)})")
    if len(top) > 1:
        notes.append(
            f"{len(top)} equally suitable RTDOSE for this plan "
            f"({', '.join(f.name for f, _ in top)}); "
            "not paired, pick one yourself"
        )
        return None
    if best < 5:  # not a plan-level physical dose
        ds = top[0][1]
        summation = str(ds.get("DoseSummationType", "?"))
        if summation.upper() in _PLAN_LEVEL:
            notes.append(
                f"the best RTDOSE is plan-level but its DoseType is "
                f"{ds.get('DoseType')}, not PHYSICAL"
            )
        else:
            notes.append(f"the best RTDOSE is not plan-level (DoseSummationType {summation})")
    return top[0][0]


def discover_with_report(patient_dir: str | Path) -> tuple[list[CaseFiles], list[str]]:
    """Cases plus the reasons any RTSTRUCT could not be placed."""
    root = Path(patient_dir)
    if not root.is_dir():
        raise FileNotFoundError(f"directory not found: {root}")

    images: dict[str, dict] = {}
    structs, plans, doses = [], [], []
    for f in _dicom_files(root):
        try:
            ds = _header(f)
        except (pydicom.errors.InvalidDicomError, OSError):
            continue
        mod = ds.get("Modality")
        if mod in IMAGE_MODALITIES:
            entry = images.setdefault(str(ds.get("SeriesInstanceUID", "")), {"files": [], "ds": ds})
            entry["files"].append(f)
        elif mod == "RTSTRUCT":
            structs.append((f, ds))
        elif mod == "RTPLAN":
            plans.append((f, ds))
        elif mod == "RTDOSE":
            doses.append((f, ds))

    structs, struct_dups = _dedupe(structs)
    plans, plan_dups = _dedupe(plans)
    doses, dose_dups = _dedupe(doses)

    dirs_of = {uid: sorted({f.parent for f in e["files"]}) for uid, e in images.items()}
    series_in_dir: dict[Path, set[str]] = {}
    for uid, ds_ in dirs_of.items():
        for d in ds_:
            series_in_dir.setdefault(d, set()).add(uid)

    cases: list[CaseFiles] = []
    skipped: list[str] = []
    for f, rs in structs:
        refs = rtstruct_referenced_series(rs)
        if not refs:
            skipped.append(
                f"{f.name}: no ReferencedFrameOfReferenceSequence, so it cannot be tied to an image"
            )
            continue
        present = [r for r in refs if r in images]
        if not present:
            skipped.append(
                f"{f.name}: its image series {[_short(r) for r in refs]} are not in this folder"
            )
            continue
        rois = [str(r.ROIName) for r in rs.get("StructureSetROISequence", [])]
        struct_uid = str(rs.SOPInstanceUID)
        plans_for = [
            (pf, pds)
            for pf, pds in plans
            if struct_uid in _ref_uids(pds, "ReferencedStructureSetSequence")
        ]
        for uid in present:
            image_ds = images[uid]["ds"]
            desc = str(image_ds.get("SeriesDescription", "") or image_ds.Modality)
            label = f"{rs.get('StructureSetLabel', '')} on {image_ds.Modality} '{desc}'"
            base_notes: list[str] = []
            if struct_uid in struct_dups:
                base_notes.append(
                    "duplicate file(s) of this structure set ignored: "
                    + ", ".join(struct_dups[struct_uid])
                )
            if len(present) > 1:
                base_notes.append(
                    f"this structure set references {len(present)} image series that are all "
                    "present; "
                    "one case per series"
                )
            image_dirs = dirs_of[uid]
            image_dir = image_dirs[0]
            if len(image_dirs) > 1:
                base_notes.append(
                    "the image series is spread over several folders; only the first is used"
                )
            if len(series_in_dir[image_dir]) > 1:
                base_notes.append(
                    "this folder holds more than one image series; the loader needs one per folder"
                )
            image_for = str(image_ds.get("FrameOfReferenceUID", ""))

            variants: list[tuple[Path | None, object | None, Path | None]] = []
            if plans_for:
                for pf, pds in plans_for:
                    plan_doses = [
                        (df, dds)
                        for df, dds in doses
                        if str(pds.SOPInstanceUID) in _ref_uids(dds, "ReferencedRTPlanSequence")
                    ]
                    dose = _choose_dose_with_notes(plan_doses, base_notes, pf, len(plans_for))
                    for uid_, names_ in (*plan_dups.items(), *dose_dups.items()):
                        if uid_ == str(pds.SOPInstanceUID) or any(
                            str(d.SOPInstanceUID) == uid_ and df_ == dose for df_, d in plan_doses
                        ):
                            base_notes.append(f"duplicate file(s) ignored: {', '.join(names_)}")
                    variants.append((pf, pds, dose))
                if len(plans_for) > 1:
                    base_notes.append(
                        f"{len(plans_for)} plans use this structure set: one case per plan"
                    )
            else:
                same_for = [
                    (df, dds)
                    for df, dds in doses
                    if str(dds.get("FrameOfReferenceUID")) == image_for
                ]
                n = len(same_for)
                reason = (
                    f"{len(plans)} RTPLAN(s) present but none references this structure set"
                    if plans
                    else "no RTPLAN"
                )
                if n:
                    nn: list[str] = []
                    dose = _choose_dose(same_for, nn)
                    base_notes.append(
                        f"{reason}; RTDOSE paired by FrameOfReferenceUID"
                        + (" (unique match)" if n == 1 else "")
                    )
                    base_notes += nn
                else:
                    dose = None
                    base_notes.append(f"{reason} and no RTDOSE shares the image's FrameOfReference")
                variants.append((None, None, dose))

            for pf, pds, dose in variants:
                notes = list(base_notes)
                if pds is not None:
                    rx, fx = _prescription(pds, notes)
                else:
                    rx, fx = None, None
                if rx is None and not any("prescription" in n or "target dose" in n for n in notes):
                    notes.append("prescription not in the data: supply prescription_gy yourself")
                if dose is not None:
                    dose_for = next(
                        (str(dds.get("FrameOfReferenceUID")) for df, dds in doses if df == dose), ""
                    )
                    if dose_for and image_for and dose_for != image_for:
                        notes.append("the RTDOSE's FrameOfReference differs from the image's")
                cases.append(
                    CaseFiles(root.name, label, image_dir, f, dose, pf, rx, fx, rois, notes)
                )
    return cases, skipped


def _choose_dose_with_notes(
    candidates, notes: list[str], plan_path: Path, n_plans: int
) -> Path | None:
    local: list[str] = []
    dose = _choose_dose(candidates, local)
    if not candidates:
        local.append(f"RTPLAN {plan_path.name} found but no RTDOSE references it")
    if n_plans > 1:
        notes.extend(f"[{plan_path.name}] {n}" for n in local)
    else:
        notes.extend(local)
    return dose


def discover_cases(patient_dir: str | Path) -> list[CaseFiles]:
    """Cases for a folder; warns about RTSTRUCTs that could not be placed (see
    ``discover_with_report`` for the reasons as data)."""
    cases, skipped = discover_with_report(patient_dir)
    for reason in skipped:
        warnings.warn(f"skipped RTSTRUCT: {reason}", stacklevel=2)
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
    elif experiment == "resolution_scan":
        cfg["resolution_scan"] = {"factors": [1, 2, 3, 4], "keys": keys}
    return cfg
