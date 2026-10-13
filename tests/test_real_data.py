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


@needs_pancreas
def test_mutual_information_recovers_known_motion_on_real_ct_anatomy_across_scanner_differences():
    """Real pancreas planning CT, 'imaged' by a different virtual scanner (coarser grid, unrelated
    frame, nonlinear intensities, noise) after a known rigid motion; start 12 mm off."""
    import glob

    from doseshift.dicomrt import load_image_series
    from doseshift.geometry import Grid, RigidMotion
    from doseshift.registration import intensity_centroid_mm, register_rigid
    from helpers import cbct_like

    g, v, _ = load_image_series(glob.glob(str(PANC / "CT_2_*"))[0])
    crop = v[:, 128:384, 128:384]
    grid = Grid(crop.shape, g.translated_voxels((128, 128, 0)).affine)
    lo, hi = grid.extent_xyz()
    centre = tuple((lo + hi) / 2)
    t, r, offset = (-4.0, 3.0, 2.0), (1.0, -1.5, 2.0), np.array([83.0, -41.0, 117.0])
    mg, mv = cbct_like(grid, np.ascontiguousarray(crop), RigidMotion(t, r, centre), offset, 1)
    t0 = intensity_centroid_mm(mg, mv) - intensity_centroid_mm(grid, crop) + [12.0, 0.0, 0.0]
    res = register_rigid(
        (grid, np.ascontiguousarray(crop)),
        (mg, mv),
        center_mm=centre,
        metric="mi",
        initial_translation_mm=tuple(t0),
    )
    assert np.allclose(res.motion.translation_mm, np.array(t) + offset, atol=1.5)
    assert np.allclose(res.motion.rotation_deg, r, atol=0.5)


@needs_pancreas
def test_real_ct_and_cbct_with_unrelated_frames_raise_instead_of_returning_identity():
    """Before the fix this returned translation 0, rotation 0 and a 1.8e308 metric."""
    import glob

    from doseshift.dicomrt import load_image_series
    from doseshift.registration import register_rigid

    plan = load_image_series(glob.glob(str(PANC / "CT_2_*"))[0])
    cbct = load_image_series(glob.glob(str(PANC / "CT_1_*"))[0])
    with pytest.raises(ValueError, match="barely overlap"):
        register_rigid(
            (plan[0], plan[1]), (cbct[0], cbct[1]), shrink_factors=(8,), smoothing_sigmas_vox=(2.0,)
        )
