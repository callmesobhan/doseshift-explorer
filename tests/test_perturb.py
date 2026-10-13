import numpy as np
import pytest

from doseshift.evaluate import evaluate
from doseshift.perturb import (
    Compose,
    DoseCoarsen,
    DoseGridOffset,
    InterpolationOrder,
    RigidShift,
    SmoothDeformation,
    StructureSampling,
)


def test_zero_perturbations_are_identity(state):
    ref = evaluate(state).flat()
    for p in (RigidShift(), DoseGridOffset((0, 0, 0)), InterpolationOrder(1), Compose()):
        out = evaluate(p(state)).flat()
        assert out == pytest.approx(ref)


def test_state_is_not_mutated(state):
    before = state.clouds["PTV"].points.copy()
    RigidShift((5, 5, 5))(state)
    DoseGridOffset((0.5, 0.5, 0.5))(state)
    assert np.array_equal(state.clouds["PTV"].points, before)
    assert state.history == ()


def test_patient_shift_equals_opposite_dose_grid_shift(state):
    """Moving the patient by +t is the same as locating the dose at -t (the core invariant)."""
    t = np.array([3.0, -6.0, 1.5])  # mm, multiples of the 3 mm dose grid -> no extra rounding
    patient = RigidShift(tuple(t))(state).structure_doses()
    voxels = tuple(-t / state.dose.grid.spacing_xyz)
    dose = DoseGridOffset(voxels)(state).structure_doses()
    for name in state.names:
        assert np.allclose(patient[name], dose[name], atol=1e-9)


def test_rotation_about_isocenter_keeps_symmetric_target_centroid(state):
    moved = RigidShift(rotation_deg=(2.0, -3.0, 5.0))(state)
    assert np.allclose(moved.clouds["PTV"].centroid, state.clouds["PTV"].centroid, atol=0.2)
    assert not np.allclose(moved.clouds["Rectum"].points, state.clouds["Rectum"].points)


def test_physical_shifts_never_improve_target_coverage(state):
    base = evaluate(state).flat()["CTV:D98%"]
    for t in ([4, 0, 0], [0, -6, 0], [0, 0, 8], [-6, 6, 0]):
        assert evaluate(RigidShift(tuple(t))(state)).flat()["CTV:D98%"] <= base + 0.15


def test_dose_cloud_metrics_ignore_setup_shift_but_not_grid_assumptions(state):
    base = evaluate(state).flat()
    shifted = evaluate(RigidShift((0, 4, 0))(state)).flat()
    assert shifted["PTV:GI"] == pytest.approx(base["PTV:GI"])
    assert shifted["PTV:CI"] < base["PTV:CI"]  # TV_PIV moves, PIV does not
    coarse = evaluate(DoseCoarsen(2)(state)).flat()
    assert coarse["PTV:GI"] != pytest.approx(base["PTV:GI"], rel=1e-3)


def test_smooth_deformation_rms_and_determinism(state):
    p = SmoothDeformation(rms_mm=3.0, corr_mm=25.0, seed=7)
    a, b = p(state), p(state)
    pts, ref = a.clouds["Rectum"].points, state.clouds["Rectum"].points
    assert np.array_equal(pts, b.clouds["Rectum"].points)
    rms = np.sqrt(((pts - ref) ** 2).sum(axis=1).mean() / 3)
    assert 0.8 < rms < 6.0
    assert np.allclose(
        SmoothDeformation(rms_mm=0.0)(state).clouds["PTV"].points, state.clouds["PTV"].points
    )


def test_dose_offset_changes_only_location(state):
    out = DoseGridOffset((0.5, 0.0, 0.0))(state)
    assert np.array_equal(out.dose.data, state.dose.data)
    assert np.isclose(out.dose.grid.origin_xyz[0] - state.dose.grid.origin_xyz[0], 1.5)


def test_coarsen_geometry_and_mean(state):
    out = DoseCoarsen(2)(state)
    assert np.allclose(out.dose.grid.spacing_xyz, 2 * state.dose.grid.spacing_xyz)
    assert out.dose.data.mean() == pytest.approx(state.dose.data.mean(), rel=1e-6)
    lo0, hi0 = state.dose.grid.extent_xyz()
    lo1, hi1 = out.dose.grid.extent_xyz()
    assert np.allclose((lo0 + hi0) / 2, (lo1 + hi1) / 2)  # still centred


def test_structure_sampling_preserves_volume_and_guards_order(state):
    v0 = state.clouds["PTV"].total_volume
    frac = StructureSampling("dose_fractional", subvoxel=3)(state).clouds["PTV"].total_volume
    assert frac == pytest.approx(v0, rel=0.03)
    with pytest.raises(ValueError):
        StructureSampling()(RigidShift((1, 0, 0))(state))


def test_history_records_order(state):
    out = Compose(DoseCoarsen(2), StructureSampling("dose_center"), RigidShift((1, 0, 0)))(state)
    assert [h.split()[0] for h in out.history] == [
        "dose-coarsen",
        "structures-on-dose-grid",
        "shift[1,0,0]mm",
    ]
