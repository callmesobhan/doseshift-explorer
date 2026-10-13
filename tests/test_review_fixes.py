"""Regression tests for the findings of the independent scientific review."""

import warnings

import numpy as np
import pytest
from scipy import ndimage
from scipy.special import ndtr

from doseshift.dvh import DVH
from doseshift.evaluate import evaluate, evaluate_doses
from doseshift.experiments import (
    _dice,
    audit,
    audit_equivalence,
    directional_sweeps,
    equivalent_shift_mm,
    margin_scan,
)
from doseshift.geometry import Grid
from doseshift.metrics import Constraint
from doseshift.montecarlo import ErrorModel, simulate
from doseshift.perturb import DoseGridOffset, SmoothDeformation
from doseshift.state import DoseGrid, MaskStructure, PlanState, PointCloud
from doseshift.synthetic import make_dose, make_state


# ---- equivalent shift: sign-aware, honest about being a lower bound ----
@pytest.fixture(scope="module")
def sweeps(state):
    return directional_sweeps(state, max_mm=10.0, step_mm=0.5)


@pytest.mark.parametrize(
    "voxels,key,direction",
    [((0, 1, 0), "Rectum:D2%", "-y"), ((0, -1, 0), "Rectum:D2%", "+y")],
)
def test_one_voxel_offset_is_recovered_at_the_right_sign(state, voxels, key, direction):
    """A 3.0 mm dose displacement is matched by a ~3 mm shift on the side that changes the metric
    the same way. Sign-agnostic matching used the steeper opposite side and reported 1.6 mm."""
    scen = {"offset": DoseGridOffset(voxels)}
    df = audit(state, scen)
    eq = audit_equivalence(state, df, [key], max_mm=10.0).iloc[0]
    assert eq["displacement_mm"] == pytest.approx(3.0)
    assert eq["direction"] == direction
    assert eq["equivalent_shift_mm"] == pytest.approx(3.0, abs=0.6)


def test_equivalent_is_a_lower_bound_for_oblique_changes_and_displacement_is_reported(state):
    df = audit(state, {"x": DoseGridOffset((1, 0, 0)), "z": DoseGridOffset((0, 0, 1))})
    eq = audit_equivalence(state, df, ["PTV:D95%"], max_mm=10.0)
    assert np.allclose(eq["displacement_mm"], 3.0)
    assert ((eq["equivalent_shift_mm"] > 0) & (eq["equivalent_shift_mm"] <= 3.3)).all()


def test_displacement_is_only_reported_for_actual_displacements(state):
    df = audit(state)
    disp = df.set_index("scenario")["displacement_mm"]
    assert np.isnan(disp["nominal"]) and np.isnan(disp["dose grid coarsened x2"])
    assert disp["centre/corner +1/2 voxel (ijk)"] == pytest.approx(0.5 * 3.0 * np.sqrt(3))


def test_no_rigid_shift_can_raise_target_coverage(sweeps):
    mm, direction = equivalent_shift_mm(sweeps, "PTV:D95%", +0.75, min_delta=0.1)
    assert np.isnan(mm) and direction == ""
    assert np.isnan(equivalent_shift_mm(sweeps, "PTV:D95%", float("nan"))[0])
    assert equivalent_shift_mm(sweeps, "PTV:D95%", -0.5)[0] > 0
    # the unsigned match still exists for comparison, and is never larger
    unsigned = equivalent_shift_mm(sweeps, "PTV:D95%", +0.75, signed=False)[0]
    assert np.isfinite(unsigned)


# ---- Paddick CI no longer depends on how the dose grid is sampled ----
def _ci_on_grid(dose_mm, sigma_p=3.2, radius=12.0):
    ct = Grid.centered((60, 60, 60), (1.0, 1.0, 1.0))
    ijk = np.indices((60, 60, 60)).reshape(3, -1).T[:, ::-1]
    mask = (np.linalg.norm(ct.index_to_world(ijk), axis=1) <= radius).reshape(60, 60, 60)
    n = int(120 / dose_mm)
    dg = Grid.centered((n, n, n), (dose_mm,) * 3)
    d = dg.index_to_world(np.indices((n, n, n)).reshape(3, -1).T[:, ::-1])
    sd = np.linalg.norm(d, axis=1) - radius  # 100 % isodose is exactly the target surface
    dose = (120.0 * ndtr(-sd / sigma_p)).reshape(n, n, n)
    st = PlanState.from_masks(DoseGrid(dg, dose), [MaskStructure("T", "ptv", mask, ct)], 60.0)
    return evaluate(st, {"T": ["CI"]}).flat()["T:CI"]


def test_conformity_index_is_nearly_independent_of_dose_grid_spacing():
    """Before the fix TV_PIV (interpolated) and PIV (node count) disagreed: CI 0.97 / 0.91 / 0.82 /
    0.58 for 1 / 2 / 3 / 5 mm grids on a plan whose isodose is exactly the target surface."""
    ci = {m: _ci_on_grid(m) for m in (1, 2, 5)}
    assert ci[1] == pytest.approx(1.0, abs=0.02)
    assert ci[2] > 0.95 and ci[5] > 0.88


# ---- blur engine ----
def test_blur_matches_the_sampler_when_sigma_is_below_one_dose_voxel(state):
    """sigma = 1 mm on the 3 mm phantom grid: filtering on the original grid under-blurred 7x."""
    model = ErrorModel((0, 0, 0), (1.0, 1.0, 1.0), n_fractions=1200)
    s = simulate(state, model, 1, seed=1, method="sampled").metrics.iloc[0]
    b = simulate(state, model, 1, seed=1, method="blur").metrics.iloc[0]
    assert b["PTV:D95%"] == pytest.approx(s["PTV:D95%"], abs=0.3)
    assert b["PTV:Dmean"] == pytest.approx(s["PTV:Dmean"], abs=0.15)


def test_blur_rejects_oblique_dose_grids(case):
    from scipy.spatial.transform import Rotation

    direction = Rotation.from_euler("z", 30, degrees=True).as_matrix()
    g = Grid.from_spacing((20, 30, 30), (3.0, 3.0, 3.0), (-45.0, -45.0, -30.0), direction)
    st = PlanState.from_masks(DoseGrid(g, np.full(g.shape, 60.0)), [case.masks["CTV"]], 60.0)
    with pytest.raises(ValueError, match="aligned with the patient axes"):
        simulate(st, ErrorModel(rand_sigma_mm=(2.0, 2.0, 2.0)), 1, method="blur")


def test_blur_maps_sigma_onto_world_axes_for_permuted_grids(case):
    """Same dose field on an axis-aligned grid and on one whose i/j axes are swapped (90 deg about
    z): blurring with anisotropic sigma must give the same result."""
    from scipy.spatial.transform import Rotation

    shape, sp = (48, 60, 60), (3.0, 3.0, 3.0)
    a = Grid.centered(shape, sp)
    rot = Rotation.from_euler("z", 90, degrees=True).as_matrix()
    centre = np.array([(shape[2] - 1) / 2, (shape[1] - 1) / 2, (shape[0] - 1) / 2])
    A = rot @ np.diag(sp)
    b = Grid.from_spacing(shape, sp, tuple(-A @ centre), rot)
    model = ErrorModel((0, 0, 0), (6.0, 2.0, 2.0), n_fractions=1)
    out = []
    for g in (a, b):
        st = PlanState.from_masks(DoseGrid(g, make_dose(g)), list(case.masks.values()), 60.0)
        out.append(simulate(st, model, 1, method="blur").metrics.iloc[0])
    for k in ("PTV:Dmean", "PTV:D95%", "Rectum:Dmean"):
        assert out[0][k] == pytest.approx(out[1][k], abs=0.35), k


# ---- Monte Carlo reporting ----
def test_probability_comes_with_a_standard_error(state):
    r = simulate(state, ErrorModel(n_fractions=2), 50, seed=0, method="blur")
    c = Constraint("CTV", "D98%", ">=", 95.0, pct_of_rx=True)
    p, se = r.probability(c), r.probability_se(c)
    assert se == pytest.approx(np.sqrt(p * (1 - p) / 50))


def test_margin_scan_reports_standard_error():
    df = margin_scan(
        lambda m: make_state(ctv_margin_mm=m),
        [5.0],
        ErrorModel((2,) * 3, (3,) * 3),
        Constraint("CTV", "Dmin", ">=", 95.0, pct_of_rx=True),
        n_courses=60,
    )
    assert {"probability", "se", "n_courses"} <= set(df.columns)
    assert df["se"].iloc[0] == pytest.approx(
        np.sqrt(df["probability"].iloc[0] * (1 - df["probability"].iloc[0]) / 60)
    )


# ---- Dice across grids ----
def _sphere(grid_spacing, r):
    g = Grid.centered(
        (int(80 / grid_spacing[2]), int(80 / grid_spacing[1]), int(80 / grid_spacing[0])),
        grid_spacing,
    )
    ijk = np.indices(g.shape).reshape(3, -1).T[:, ::-1]
    m = (np.linalg.norm(g.index_to_world(ijk), axis=1) <= r).reshape(g.shape)
    return MaskStructure("S", "oar", m, g)


@pytest.mark.parametrize("r", [3.0, 6.0, 12.0])
def test_dice_is_symmetric_and_bounded_across_different_grids(r):
    a, b = _sphere((2.0, 1.0, 1.0), r), _sphere((3.0, 1.5, 1.5), r)
    ab, ba = _dice(a, b), _dice(b, a)
    assert ab == pytest.approx(ba, abs=1e-12)  # one-sided estimates gave 0.72 vs 1.10 at r = 3
    assert 0.0 < ab <= 1.0
    assert ab > (0.6 if r == 3.0 else 0.8)


# ---- edge cases ----
def test_nan_dose_is_rejected_not_averaged_into_plausible_numbers():
    with pytest.raises(ValueError, match="finite"):
        DVH(np.array([10.0, float("nan"), 30.0]), np.ones(3))


def test_dose_at_cc_is_nan_when_the_structure_is_smaller_than_requested():
    d = DVH(np.array([10.0, 20.0]), np.array([500.0, 500.0]))  # 1 cc in total
    assert np.isnan(d.dose_at_volume_cc(5.0)) and np.isfinite(d.dose_at_volume_cc(0.5))


def test_geud_a_zero_is_the_geometric_mean():
    d = DVH(np.array([10.0, 40.0]), np.array([1.0, 1.0]))
    assert d.generalized_eud(0.0) == pytest.approx(20.0)
    assert d.generalized_eud(1e-6) == pytest.approx(20.0, rel=1e-4)


def test_structures_with_very_few_sample_points_warn(state):
    pts = state.clouds["Rectum"].points[:5]
    tiny = state.evolve(clouds={**state.clouds, "Rectum": PointCloud(pts, np.full(5, 12.0))})
    doses = tiny.structure_doses()
    with pytest.warns(UserWarning, match="only 5 sample points"):
        evaluate_doses(tiny, doses, {"Rectum": ["Dmean"]})


@pytest.mark.parametrize("seed", [0, 3, 11])
def test_smooth_deformation_realises_the_requested_rms_and_has_no_rigid_part(state, seed):
    out = SmoothDeformation(rms_mm=2.0, corr_mm=25.0, seed=seed)(state)
    ref = np.concatenate([state.clouds[n].points for n in state.names])
    new = np.concatenate([out.clouds[n].points for n in state.names])
    disp = new - ref
    assert np.sqrt((disp**2).mean()) == pytest.approx(
        2.0, rel=0.02
    )  # was 0.65-1.25 mm per component
    assert np.abs(disp.mean(axis=0)).max() < 1e-9


def test_structure_sampling_scenarios_keep_structure_volume(state):
    from doseshift.perturb import StructureSampling

    for p in (StructureSampling("dose_center"), StructureSampling("dose_fractional")):
        for n in state.names:
            v0, v1 = state.clouds[n].total_volume, p(state).clouds[n].total_volume
            assert abs(v1 - v0) / v0 < 0.08, (p.mode, n)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # a well-covered phantom raises no volume warnings
        ndimage.label(np.ones((2, 2)))
