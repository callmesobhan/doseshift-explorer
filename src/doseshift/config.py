"""YAML experiment configs -> library objects. See experiments/*.yaml for annotated examples."""

from __future__ import annotations

from pathlib import Path

import yaml

from .geometry import Grid
from .metrics import Constraint
from .montecarlo import ErrorModel
from .state import PlanState

EXPERIMENTS = (
    "sweep",
    "grid_audit",
    "monte_carlo",
    "margin_scan",
    "contour_variation",
    "resolution_scan",
)


def load_config(path: str | Path) -> dict:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"config not found: {path}")
    cfg = yaml.safe_load(path.read_text())
    if isinstance(cfg, dict) and "config" in cfg and "doseshift" in cfg:
        cfg = cfg["config"]  # a provenance.json written by a previous run re-runs it
    validate_config(cfg)
    return cfg


def validate_config(cfg: dict) -> None:
    if not isinstance(cfg, dict):
        raise ValueError("config must be a mapping")
    if cfg.get("experiment") not in EXPERIMENTS:
        raise ValueError(f"'experiment' must be one of {EXPERIMENTS}")
    case = cfg.get("case")
    if not isinstance(case, dict) or case.get("type") not in ("synthetic", "dicom"):
        raise ValueError("case.type must be 'synthetic' or 'dicom'")
    if case["type"] == "dicom":
        missing = {"dose", "rtstruct", "prescription_gy", "roles"} - set(case)
        if "ct_dir" not in case and "image_dir" not in case:
            missing.add("image_dir")
        if missing:
            raise ValueError(f"dicom case is missing {sorted(missing)}")
    mc = cfg.get("monte_carlo", {})
    if mc and int(mc.get("n_courses", 1)) < 1:
        raise ValueError("monte_carlo.n_courses must be >= 1")
    if cfg["experiment"] == "contour_variation":
        if case["type"] != "dicom":
            raise ValueError("contour_variation needs real contours (case.type: dicom)")
        cv = cfg.get("contour_variation", {})
        for alt in cv.get("alternatives", []) or [None]:
            if (
                not alt
                or not {"label", "rtstruct", "names"} <= set(alt)
                or not ("image_dir" in alt or "ct_dir" in alt)
            ):
                raise ValueError(
                    "each contour_variation alternative needs label, image_dir, "
                    "rtstruct and names (planning ROI -> alternative ROI)"
                )
    if cfg["experiment"] == "margin_scan" and case["type"] != "synthetic":
        raise ValueError("margin_scan needs the synthetic phantom (it re-draws the CTV)")


def image_dir(case: dict) -> str:
    """The CT/MR series directory (``image_dir``; ``ct_dir`` is accepted as an alias)."""
    return case.get("image_dir", case.get("ct_dir"))


def build_state(case: dict) -> tuple[PlanState, tuple[Grid, object] | None]:
    """Returns the PlanState and an optional (grid, HU) CT for display."""
    if case["type"] == "synthetic":
        from .synthetic import make_case

        c = make_case(int(case.get("seed", 0)), float(case.get("ctv_margin_mm", 5.0)))
        return c.to_state(), (c.ct_grid, c.ct_hu)
    from .dicomrt import load_plan_state_with_image

    for key, path in (
        ("image_dir", image_dir(case)),
        ("dose", case["dose"]),
        ("rtstruct", case["rtstruct"]),
    ):
        if not Path(path).exists():
            raise FileNotFoundError(
                f"{key}: {path} not found (relative paths resolve against the current directory; "
                "run from the repository root, and fetch data with scripts/download_tcia.py)"
            )
    state, image = load_plan_state_with_image(
        image_dir(case),
        case["dose"],
        case["rtstruct"],
        case["roles"],
        case["prescription_gy"],
        "reference dose" if case.get("prescription_is_reference") else "Rx",
    )
    return state, image


def error_model_from(d: dict | None) -> ErrorModel:
    d = dict(d or {})
    for k in ("sys_sigma_mm", "rand_sigma_mm", "sys_sigma_deg", "rand_sigma_deg"):
        if k in d:
            d[k] = tuple(float(v) for v in d[k])
    return ErrorModel(**d)


def constraints_from(items: list[dict] | None) -> list[Constraint]:
    return [Constraint(**c) for c in (items or [])]
