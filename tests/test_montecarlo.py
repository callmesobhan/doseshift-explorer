import numpy as np
import pytest

from doseshift.metrics import Constraint
from doseshift.montecarlo import ErrorModel, simulate

ZERO = ErrorModel((0, 0, 0), (0, 0, 0), n_fractions=3)


def test_van_herk_recipe_numbers():
    m = ErrorModel((2.0, 2.0, 2.0), (3.0, 3.0, 3.0))
    assert m.van_herk_margin_mm() == pytest.approx(6.9456, abs=1e-3)
    assert ErrorModel((0, 0, 0), (0, 0, 0)).van_herk_margin_mm().max() == 0.0


@pytest.mark.parametrize("method", ["sampled", "blur"])
def test_zero_error_reproduces_nominal(state, method):
    r = simulate(state, ZERO, n_courses=3, seed=0, method=method)
    nominal = r.nominal.flat()
    for key in ("PTV:D95%", "CTV:D98%", "Rectum:Dmean"):
        assert r.metrics[key].to_numpy() == pytest.approx(nominal[key], abs=1e-6)


def test_seeded_runs_are_reproducible(state):
    m = ErrorModel(n_fractions=4)
    a = simulate(state, m, n_courses=5, seed=11).metrics
    b = simulate(state, m, n_courses=5, seed=11).metrics
    c = simulate(state, m, n_courses=5, seed=12).metrics
    assert a.equals(b) and not a.equals(c)


def test_blur_engine_matches_sampler_in_many_fraction_limit(state):
    model = ErrorModel((0, 0, 0), (3.0, 3.0, 3.0), n_fractions=1500)
    s = simulate(state, model, n_courses=1, seed=3, method="sampled").metrics.iloc[0]
    b = simulate(state, model, n_courses=1, seed=3, method="blur").metrics.iloc[0]
    assert s["PTV:Dmean"] == pytest.approx(b["PTV:Dmean"], abs=0.3)
    assert s["PTV:D95%"] == pytest.approx(b["PTV:D95%"], abs=0.8)
    assert s["Rectum:Dmean"] == pytest.approx(b["Rectum:Dmean"], abs=0.2)


def test_blur_rejects_random_rotations(state):
    with pytest.raises(ValueError):
        simulate(state, ErrorModel(rand_sigma_deg=(1, 1, 1)), n_courses=1, method="blur")


def test_larger_systematic_error_lowers_probability_of_success(state):
    c = Constraint("CTV", "D98%", ">=", 95.0, pct_of_rx=True)
    p = [
        simulate(
            state, ErrorModel((s,) * 3, (2.0,) * 3, n_fractions=5), 120, seed=5, method="blur"
        ).probability(c)
        for s in (0.5, 3.0, 6.0)
    ]
    assert p[0] >= p[1] >= p[2] and p[0] > 0.95 and p[2] < 0.8


def test_summary_and_band_shapes(state):
    r = simulate(state, ErrorModel(n_fractions=3), n_courses=6, seed=0)
    assert r.band("PTV").shape == (3, len(r.dvh_axis))
    assert {"nominal", "mean", "p5", "p95"} <= set(r.summary().columns)
    assert np.all(r.dvh_curves["PTV"][:, 0] == 1.0)


def test_reproduces_van_herk_margin_recipe():
    """The phantom matches the recipe's assumptions (95 % isodose on the PTV, sigma_p = 3.2 mm),
    so P(CTV Dmin >= 95 % Rx) must be ~0.9 at the recipe margin and fall steeply below it."""
    from doseshift.experiments import margin_scan

    c = Constraint("CTV", "Dmin", ">=", 95.0, pct_of_rx=True)
    full = ErrorModel((2.0,) * 3, (3.0,) * 3)
    p3 = margin_scan(make_state_at, [3.0, full.van_herk_margin_mm()[0]], full, c, 800, 4)
    assert p3.loc[0, "probability"] < 0.1
    assert 0.84 <= p3.loc[1, "probability"] <= 0.95  # recipe target 0.90

    one_d = ErrorModel((2.0, 0, 0), (3.0, 0, 0))
    m1 = one_d.van_herk_margin_mm(systematic_factor=1.2816)[0]  # ~4.5 mm
    p1 = margin_scan(make_state_at, [m1], one_d, c, 800, 4)
    assert 0.84 <= p1.loc[0, "probability"] <= 0.96


def make_state_at(margin: float):
    from doseshift.synthetic import make_state

    return make_state(ctv_margin_mm=margin)
