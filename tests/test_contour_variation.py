import numpy as np
import pytest
from scipy import ndimage

from doseshift.experiments import contour_variation
from doseshift.state import MaskStructure


def _alt(case, name, mask):
    s = case.masks[name]
    return {name: MaskStructure(name, s.role, mask, s.grid)}


def test_a_known_contour_shift_is_recovered_in_mm(case, state):
    """Rolling the rectum mask 2 voxels (4 mm) along -y is the planning contour moved 4 mm."""
    rect = case.masks["Rectum"]
    moved = np.roll(rect.mask, -2, axis=1)  # axis 1 = y; -y brings it towards the PTV
    df = contour_variation(
        state, {"shifted": _alt(case, "Rectum", moved)}, metrics={"Rectum": ["Dmean", "D2%"]}
    )
    row = df[df["metric"] == "Dmean"].iloc[0]
    assert row["rigid_comparable"] and row["volume_ratio"] == pytest.approx(1.0)
    assert row["centroid_shift_mm"] == pytest.approx(4.0, abs=0.1)
    assert row["centroid_dy_mm"] == pytest.approx(-4.0, abs=0.1)
    assert 0.3 < row["dice"] < 0.95
    assert row["delta"] > 0  # dose rises towards the target
    assert row["equivalent_shift_mm"] == pytest.approx(4.0, abs=1.0)


def test_identical_contours_change_nothing(case, state):
    df = contour_variation(state, {"same": _alt(case, "Rectum", case.masks["Rectum"].mask)})
    assert (df["delta"].abs() < 1e-9).all() and (df["dice"] == 1.0).all()
    assert (df["equivalent_shift_mm"] == 0).all()


def test_changed_organ_volume_is_flagged_not_rigid_comparable(case, state):
    grown = ndimage.binary_dilation(case.masks["Rectum"].mask, iterations=3)
    df = contour_variation(state, {"grown": _alt(case, "Rectum", grown)})
    assert not df["rigid_comparable"].any() and df["volume_ratio"].iloc[0] > 1.25
