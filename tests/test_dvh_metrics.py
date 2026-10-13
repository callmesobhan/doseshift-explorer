import numpy as np
import pytest

from doseshift.dvh import DVH
from doseshift.metrics import Constraint, MetricContext, compute_metric, metric_unit


def ctx(doses, volumes, rx=60.0, piv=0.0, v50=0.0):
    d = DVH(doses, volumes)
    return MetricContext(d, np.asarray(doses, float), np.asarray(volumes, float), rx, piv, v50)


def test_constant_dose():
    d = DVH(np.full(10, 42.0), np.ones(10))
    assert d.dose_at_volume_fraction(0.95) == 42.0
    assert d.dmean == d.dmin == d.dmax == 42.0
    assert d.volume_fraction_at_dose(42.0) == 1.0
    assert d.volume_fraction_at_dose(42.1) == 0.0


def test_weighted_dose_at_volume():
    d = DVH([10.0, 20.0], [1.0, 3.0])  # 75 % of the volume receives 20 Gy
    assert d.dose_at_volume_fraction(0.5) == 20.0
    assert d.dose_at_volume_fraction(0.75) == 20.0
    assert d.dose_at_volume_fraction(0.9) == 10.0
    assert d.dmean == pytest.approx(17.5)


def test_extremes_and_monotone_curve():
    rng = np.random.default_rng(1)
    d = DVH(rng.uniform(0, 60, 500), rng.uniform(1, 2, 500))
    assert d.dose_at_volume_fraction(0.0) == d.dmax
    assert d.dose_at_volume_fraction(1.0) == d.dmin
    c = d.curve(np.linspace(0, 70, 200))
    assert c[0] == 1.0 and c[-1] == 0.0 and np.all(np.diff(c) <= 1e-12)


def test_dose_at_cc_matches_fraction():
    d = DVH(np.arange(100, dtype=float), np.full(100, 100.0))  # 10 cc total
    assert d.total_volume_cc == pytest.approx(10.0)
    assert d.dose_at_volume_cc(1.0) == d.dose_at_volume_fraction(0.1)


def test_geud_a1_is_mean_and_increases_towards_max_with_a():
    rng = np.random.default_rng(2)
    d = DVH(rng.uniform(5, 50, 300), np.ones(300))
    assert d.generalized_eud(1.0) == pytest.approx(d.dmean)
    values = [d.generalized_eud(a) for a in (1, 4, 16, 60)]
    assert values == sorted(values) and values[-1] < d.dmax
    assert d.generalized_eud(-8.0) < d.dmean  # negative a is cold-spot weighted


@pytest.mark.parametrize(
    "spec,expected",
    [
        ("Dmean", 30.0),
        ("D50%", 50.0),
        ("D51%", 10.0),
        ("V30Gy", 50.0),
        ("V50%", 50.0),
        ("V30Gy_cc", 1.0),
    ],
)
def test_metric_grammar(spec, expected):
    c = ctx([10.0, 50.0], [1000.0, 1000.0])
    assert compute_metric(spec, c) == pytest.approx(expected)


def test_unknown_metric_raises():
    with pytest.raises(ValueError):
        compute_metric("Q95", ctx([1.0], [1.0]))


def test_units():
    assert [metric_unit(s) for s in ("D95%", "V20Gy", "V30Gy_cc", "HI", "gEUD(a=8)")] == [
        "Gy",
        "%",
        "cc",
        "-",
        "Gy",
    ]


def test_paddick_ci_perfect_and_partial():
    # target exactly equals the prescription isodose -> CI = 1
    c = ctx([60.0] * 4, [1.0] * 4, piv=4.0)
    assert compute_metric("CI", c) == pytest.approx(1.0)
    # half of the target is cold, isodose volume unchanged -> CI = 0.5^2 / ... = 0.25 * 4/(4) ...
    c = ctx([60.0, 60.0, 30.0, 30.0], [1.0] * 4, piv=4.0)
    assert compute_metric("CI", c) == pytest.approx(2.0**2 / (4.0 * 4.0))


def test_constraint_passes():
    c = Constraint("CTV", "D98%", ">=", 95.0, pct_of_rx=True)
    assert list(c.passes(np.array([56.0, 57.0, 58.0]), 60.0)) == [False, True, True]
    assert c.threshold(60.0) == pytest.approx(57.0)
