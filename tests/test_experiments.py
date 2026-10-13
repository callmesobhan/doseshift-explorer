import numpy as np
import pytest

from doseshift.experiments import (
    audit,
    audit_equivalence,
    directional_sweeps,
    equivalent_shift_mm,
    translation_sweep,
)
from doseshift.perturb import DoseGridOffset


def test_translation_sweep_zero_row_is_nominal(state):
    df = translation_sweep(state, "y", [0.0, 3.0])
    assert df.loc[0, "shift_mm"] == 0.0
    assert df.loc[1, "PTV:D95%"] < df.loc[0, "PTV:D95%"]


def test_audit_layout(state):
    df = audit(state, {"half": DoseGridOffset((0.5, 0.5, 0.5))})
    assert list(df["scenario"]) == ["nominal", "half"]
    assert list(df["kind"]) == ["-", "numerical"]


def test_equivalent_shift_is_monotone_in_delta(state):
    sweeps = directional_sweeps(state, max_mm=10.0, step_mm=0.5)
    small, _ = equivalent_shift_mm(sweeps, "PTV:D95%", -0.5)
    big, _ = equivalent_shift_mm(sweeps, "PTV:D95%", -5.0)
    assert 0 < small < big < np.inf
    assert equivalent_shift_mm(sweeps, "PTV:D95%", 0.0) == (0.0, "")
    assert equivalent_shift_mm(sweeps, "PTV:D95%", -1e6)[0] == np.inf


def test_equivalence_roundtrip_recovers_a_known_shift(state):
    """A pure grid offset of d mm along y should be matched by about d mm of patient shift."""
    d_mm = 3.0
    scenario = {"dose +1 voxel (y)": DoseGridOffset((0.0, d_mm / 3.0, 0.0))}
    a = audit(state, scenario)
    eq = audit_equivalence(state, a, ["PTV:D95%"])
    assert eq.loc[0, "equivalent_shift_mm"] == pytest.approx(d_mm, abs=0.6)
