"""Checks against real public data (TCIA, CC BY 4.0). Skipped when the files are not downloaded:
see data/README.md and scripts/download_tcia.py."""

from pathlib import Path

import numpy as np
import pytest

from doseshift.dicomrt import discover_cases, load_plan_state
from doseshift.evaluate import evaluate

RAW = Path(__file__).parent.parent / "data" / "raw"
PANC = RAW / "pancreas-ct-cbct-seg" / "Pancreas-CT-CB_014"
VS = RAW / "vs-seg" / "VS-SEG-135"
needs_pancreas = pytest.mark.skipif(not PANC.exists(), reason="pancreas data not downloaded")
needs_vs = pytest.mark.skipif(not VS.exists(), reason="VS-SEG data not downloaded")


@needs_pancreas
def test_pancreas_discovery_pairs_planning_ct_with_planning_contours_and_dose():
    cases = discover_cases(PANC)
    plan = next(c for c in cases if "planCT" in " ".join(c.roi_names))
    assert "DI" in plan.label and plan.rtdose is not None and plan.prescription_gy is None
    assert len(cases) == 3  # planning + two CBCT-day structure sets


@needs_pancreas
def test_pancreas_real_dose_grid_differs_from_ct_grid_and_loads():
    plan = next(c for c in discover_cases(PANC) if "planCT" in " ".join(c.roi_names))
    state = load_plan_state(
        plan.image_dir,
        plan.rtdose,
        plan.rtstruct,
        {"ROI": "ptv", "Stomach_duo_planCT": "oar"},
        46.0,
        rx_label="reference dose",
    )
    assert state.dose.grid.shape != state.masks["ROI"].grid.shape
    assert np.allclose(state.dose.grid.spacing_xyz, (2.5, 2.5, 3.0))
    d2 = evaluate(state, {"Stomach_duo_planCT": ["D2%"], "ROI": ["D95%"]}).flat()
    assert d2["Stomach_duo_planCT:D2%"] == pytest.approx(48.6, abs=1.0)


@needs_vs
def test_vs_seg_rtplan_prescription_and_negative_frame_offsets():
    (case,) = discover_cases(VS)
    assert case.prescription_gy == 13.0 and case.n_fractions == 1
    state = load_plan_state(case.image_dir, case.rtdose, case.rtstruct, {"TV": "ptv"}, 13.0)
    assert state.dose.grid.affine[2, 2] < 0 or state.dose.grid.affine[:3, 2].sum() < 0
    ev = evaluate(state, {"TV": ["D95%"]}).flat()
    assert ev["TV:D95%"] == pytest.approx(13.3, abs=0.7)  # plan covers the tumour at ~Rx
    assert state.clouds["TV"].total_volume / 1000 == pytest.approx(1.0, abs=0.3)  # cc
