import numpy as np
import pytest

from doseshift.dicomrt import (
    alignment_report,
    load_ct,
    load_plan_state,
    load_rtdose,
    load_rtstruct_masks,
    write_case,
)
from doseshift.evaluate import evaluate

ROLES = {"PTV": "ptv", "CTV": "ctv", "Rectum": "oar", "Bladder": "oar"}


@pytest.fixture(scope="module")
def written(tmp_path_factory, case):
    root = tmp_path_factory.mktemp("dcm")
    return write_case(case, root), write_case(case, root / "rev", reverse_dose_z=True)


def test_ct_roundtrip(case, written):
    grid, hu, _ = load_ct(written[0]["ct_dir"])
    assert grid.allclose(case.ct_grid, atol=1e-4)
    assert np.abs(hu - case.ct_hu).max() <= 0.5 + 1e-3  # integer HU quantisation


@pytest.mark.parametrize("which", [0, 1], ids=["ascending", "reversed-frames"])
def test_dose_roundtrip_including_reversed_frame_order(case, written, which):
    dose, _ = load_rtdose(written[which]["dose"])
    probe = case.ct_grid.mask_points(case.masks["PTV"].mask)[::97]
    ref = case.dose
    from doseshift.geometry import sample_volume

    a = sample_volume(ref.grid, ref.data, probe)
    b = sample_volume(dose.grid, dose.data, probe)
    assert np.allclose(a, b, atol=1e-4)
    if which == 0:
        assert dose.grid.allclose(ref.grid, atol=1e-4)


def test_rtstruct_roundtrip_is_exact(case, written):
    masks, _ = load_rtstruct_masks(written[0]["rtstruct"], case.ct_grid)
    for name, s in case.masks.items():
        assert np.array_equal(masks[name], s.mask), name


def test_full_state_matches_in_memory_metrics(case, written):
    paths = written[0]
    state = load_plan_state(paths["ct_dir"], paths["dose"], paths["rtstruct"], ROLES, 60.0)
    ref = evaluate(case.to_state()).flat()
    got = evaluate(state).flat()
    for key in ("PTV:D95%", "CTV:D98%", "Rectum:Dmean", "PTV:CI", "PTV:GI"):
        assert got[key] == pytest.approx(ref[key], rel=1e-4, abs=1e-4)


def test_missing_roi_is_reported(written):
    paths = written[0]
    with pytest.raises(KeyError):
        load_plan_state(paths["ct_dir"], paths["dose"], paths["rtstruct"], {"Nope": "oar"}, 60.0)


def test_alignment_report_mentions_subvoxel_phase(written):
    paths = written[0]
    text = "\n".join(alignment_report(paths["ct_dir"], paths["dose"], paths["rtstruct"]))
    assert "CT == DOSE : True" in text and "sub-voxel" in text
