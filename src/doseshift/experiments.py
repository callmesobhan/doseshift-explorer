"""Deterministic experiments: translation sweeps and the grid-assumption audit.

The audit's headline output is ``equivalent_shift_mm``: how large a *physical* setup error
produces the same change in a metric as a *numerical* convention choice. It puts "we assumed the
wrong voxel convention" and "the patient was misaligned" on one axis.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import pandas as pd

from .evaluate import evaluate
from .metrics import Constraint, metric_unit
from .montecarlo import ErrorModel, simulate
from .perturb import (
    Compose,
    DoseCoarsen,
    DoseGridOffset,
    InterpolationOrder,
    Perturbation,
    RigidShift,
    StructureSampling,
)
from .perturb.grid import mask_contains
from .state import MaskStructure, PlanState

AXES = {"x": 0, "y": 1, "z": 2}


def translation_sweep(
    state: PlanState,
    axis: str,
    values_mm: Sequence[float],
    metrics: dict[str, list[str]] | None = None,
) -> pd.DataFrame:
    """One row per shift of the patient along ``axis``; columns are ``STRUCT:metric``."""
    rows = []
    for v in values_mm:
        t = [0.0, 0.0, 0.0]
        t[AXES[axis]] = float(v)
        flat = evaluate(RigidShift(tuple(t))(state), metrics).flat()
        rows.append({"axis": axis, "shift_mm": float(v), **flat})
    return pd.DataFrame(rows)


def default_audit_scenarios() -> dict[str, Perturbation]:
    """Plausible grid-handling differences, each applied to an otherwise perfect alignment."""
    half = 0.5
    sc: dict[str, Perturbation] = {}
    for sign, tag in ((+1, "+"), (-1, "-")):
        sc[f"centre/corner {tag}1/2 voxel (ijk)"] = DoseGridOffset((sign * half,) * 3)
        for ax, name in enumerate("ijk"):  # index axes, not world axes (see DoseGridOffset)
            v = [0.0, 0.0, 0.0]
            v[ax] = sign * half
            sc[f"centre/corner {tag}1/2 voxel ({name})"] = DoseGridOffset(tuple(v))
    sc["dose grid coarsened x2"] = DoseCoarsen(2)
    sc["interpolation: nearest"] = InterpolationOrder(0)
    sc["interpolation: cubic"] = InterpolationOrder(3)
    sc["structures on dose grid (centre)"] = StructureSampling("dose_center")
    sc["structures on dose grid (fractional)"] = StructureSampling("dose_fractional")
    sc["dose coarsened x2 + centre sampling"] = Compose(
        DoseCoarsen(2), StructureSampling("dose_center")
    )
    return sc


def audit(
    state: PlanState,
    scenarios: dict[str, Perturbation] | None = None,
    metrics: dict[str, list[str]] | None = None,
) -> pd.DataFrame:
    """Rows = scenarios (plus ``nominal``), columns = ``STRUCT:metric``, values = metric value."""
    scenarios = default_audit_scenarios() if scenarios is None else scenarios
    rows = [
        {
            "scenario": "nominal",
            "kind": "-",
            "displacement_mm": float("nan"),
            **evaluate(state, metrics).flat(),
        }
    ]
    for label, p in scenarios.items():
        rows.append(
            {
                "scenario": label,
                "kind": p.kind,
                "displacement_mm": grid_displacement_mm(state, p),
                **evaluate(p(state), metrics).flat(),
            }
        )
    return pd.DataFrame(rows)


def grid_displacement_mm(state: PlanState, p: Perturbation) -> float:
    """Size of the physical displacement a ``DoseGridOffset`` amounts to, ``|A @ v|`` in mm.

    NaN for every other scenario: coarsening, interpolation and structure sampling are not
    displacements, which is exactly where an equivalent rigid shift adds information.
    """
    if isinstance(p, DoseGridOffset):
        return float(np.linalg.norm(state.dose.grid.affine[:3, :3] @ np.asarray(p.voxels, float)))
    return float("nan")


Sweeps = dict[str, pd.DataFrame]  # "+x", "-x", ... -> translation_sweep table


def directional_sweeps(
    state: PlanState,
    max_mm: float = 12.0,
    step_mm: float = 0.25,
    metrics: dict[str, list[str]] | None = None,
) -> Sweeps:
    """Translation sweeps in all six axis directions, reusable across many equivalence lookups."""
    values = np.arange(0.0, max_mm + step_mm, step_mm)
    return {
        f"{'+' if sign > 0 else '-'}{axis}": translation_sweep(state, axis, sign * values, metrics)
        for axis in "xyz"
        for sign in (+1, -1)
    }


def equivalent_shift_mm(
    sweeps: Sweeps, key: str, delta: float, min_delta: float = 1e-9, signed: bool = True
) -> tuple[float, str]:
    """Smallest axis-aligned rigid shift that changes ``key`` by ``delta``.

    This is a **lower bound** on the displacement behind a change, not an estimate of it: the
    metric responds more steeply to a shift in some directions than others, and only the six
    axis-aligned directions are searched. A change that really comes from displacing the dose
    by d mm along an oblique direction is reported as <= d.

    ``signed=True`` (default) only accepts shifts that move the metric the *same way* as
    ``delta``; a numerical assumption that raises target coverage then has no equivalent
    (NaN), because no rigid shift raises it. ``signed=False`` matches ``|delta|`` in either
    direction (the shallow-sided match can come from the opposite side and understate the shift).

    Returns ``(mm, direction)``: ``(0, "")`` when ``|delta| < min_delta`` (below the noise
    floor, meaning negligible); ``(inf, "")`` if the right-signed change is reachable but not
    within the swept range; ``(nan, "")`` if ``delta`` is not finite or no shift can change the
    metric in that direction by at least ``min_delta``.
    """
    if not np.isfinite(delta):
        return float("nan"), ""
    target = abs(delta)
    if target < max(min_delta, 1e-9):
        return 0.0, ""
    sign = np.sign(delta) if signed else None
    best, best_dir, reachable = np.inf, "", False
    for direction, df in sweeps.items():
        mm = df["shift_mm"].abs().to_numpy()
        change = df[key].to_numpy() - df[key].iloc[0]
        change = np.abs(change) if sign is None else change * sign
        reachable |= bool(np.nanmax(change) >= max(min_delta, 1e-6))  # a shift can go that way
        hit = np.nonzero(change >= target)[0]
        if hit.size == 0:
            continue
        i = hit[0]
        x = (
            mm[i]
            if i == 0
            else float(np.interp(target, [change[i - 1], change[i]], [mm[i - 1], mm[i]]))
        )
        if x < best:
            best, best_dir = x, direction
    if not reachable:
        return float("nan"), ""
    return float(best), best_dir


def noise_floor(metric: str, min_delta: float) -> float:
    """``min_delta`` is in Gy / cc / % for those metrics; dimensionless indices (HI, CI, GI) get a
    twentieth of it, since a change of 0.1 is enormous for a number near 1."""
    return min_delta * 0.05 if metric_unit(metric) == "-" else min_delta


def audit_equivalence(
    state: PlanState,
    audit_df: pd.DataFrame,
    keys: Sequence[str],
    metrics: dict[str, list[str]] | None = None,
    max_mm: float = 12.0,
    min_delta: float = 0.1,
) -> pd.DataFrame:
    """Long table: scenario x metric -> metric change and the equivalent rigid shift in mm.

    ``equivalent_shift_mm`` is a sign-aware lower bound (see ``equivalent_shift_mm``). For
    ``DoseGridOffset`` scenarios ``displacement_mm`` holds the true displacement ``|A v|``; there
    the
    equivalent is a consistency check on the sweep, not a finding. For coarsening, interpolation
    and structure sampling there is no displacement; the equivalent is the informative number.

    Changes smaller than the noise floor (``min_delta`` = 0.1 Gy / cc / % by default, one twentieth
    of that for dimensionless indices) are reported as 0 mm, meaning "negligible", not "zero".
    """
    sweeps = directional_sweeps(state, max_mm=max_mm, metrics=metrics)
    nominal = audit_df.iloc[0]
    rows = []
    for _, r in audit_df.iloc[1:].iterrows():
        for key in keys:
            delta = float(r[key] - nominal[key])
            mm, direction = equivalent_shift_mm(
                sweeps, key, delta, noise_floor(key.split(":", 1)[1], min_delta)
            )
            rows.append(
                {
                    "scenario": r["scenario"],
                    "kind": r["kind"],
                    "metric": key,
                    "delta": delta,
                    "equivalent_shift_mm": mm,
                    "direction": direction,
                    "displacement_mm": float(r.get("displacement_mm", float("nan"))),
                }
            )
    return pd.DataFrame(rows)


def margin_scan(
    state_factory: Callable[[float], PlanState],
    margins_mm: Sequence[float],
    model: ErrorModel,
    constraint: Constraint,
    n_courses: int = 500,
    seed: int = 0,
    method: str = "blur",
) -> pd.DataFrame:
    """P(constraint met) versus CTV->PTV margin; ``state_factory(margin)`` builds each plan.

    With the synthetic phantom this lands near van Herk's 90 % target at the recipe margin
    (P(CTV Dmin >= 95 % Rx) ~ 0.89 at ``model.van_herk_margin_mm()`` in 3D), but that number moves
    between about 0.79 and 0.91 with mask/dose resolution, the plateau level and the fraction model,
    so read it as consistency with the recipe, not a validation. ``method="blur"`` is the
    infinite-fraction limit and ignores ``model.n_fractions``; with 20 explicit fractions
    (``method="sampled"``) P is lower. ``se`` is the binomial standard error over courses.
    """
    metrics = {constraint.structure: [constraint.metric]}
    rows = []
    for m in margins_mm:
        r = simulate(state_factory(float(m)), model, n_courses, seed, method, metrics)
        rows.append(
            {
                "margin_mm": float(m),
                "probability": r.probability(constraint),
                "se": r.probability_se(constraint),
                "n_courses": n_courses,
            }
        )
    return pd.DataFrame(rows)


def _dice(a: MaskStructure, b: MaskStructure) -> float:
    """Symmetric Dice of two masks that may live on different grids.

    The overlap volume is estimated from both sides (A's voxels found in B, B's voxels found in A)
    and averaged, so ``dice(a, b) == dice(b, a)``; one-sided estimates differ by several percent
    (and can exceed 1) when the grids differ.
    """
    i_ab = mask_contains(b, a.grid.mask_points(a.mask)).sum() * a.grid.voxel_volume
    i_ba = mask_contains(a, b.grid.mask_points(b.mask)).sum() * b.grid.voxel_volume
    va = a.mask.sum() * a.grid.voxel_volume
    vb = b.mask.sum() * b.grid.voxel_volume
    return float(min(1.0, (i_ab + i_ba) / (va + vb))) if va + vb else float("nan")


def contour_variation(
    state: PlanState,
    alternatives: dict[str, dict[str, MaskStructure]],
    metrics: dict[str, list[str]] | None = None,
    max_mm: float = 12.0,
    min_delta: float = 0.1,
    max_volume_change: float = 0.25,
) -> pd.DataFrame:
    """Measured geometry change: the planned dose against alternative contours of the same organs.

    ``alternatives[label][name]`` is a MaskStructure named like the planning structure ``name``
    (e.g. day-of-treatment CBCT contours in the planning frame of reference). Dose stays fixed, so
    every metric change is purely geometric; each is also expressed as the rigid shift that would
    change the metric as much (``equivalent_shift_mm``, from sweeps of the planning contours).
    One row per alternative x structure x metric, with overlap descriptors (Dice, volumes,
    centroid displacement).

    ``rigid_comparable`` is False when the organ volume changed by more than
    ``max_volume_change`` (fractional): the contour *definition or extent* then differs, not just
    the position, and a rigid-equivalent shift is not a meaningful summary of that row.
    """
    names = list(next(iter(alternatives.values())))
    sub = PlanState.from_masks(state.dose, [state.masks[n] for n in names], state.prescription_gy)
    wanted = {n: (metrics or {}).get(n) for n in names}
    wanted = {n: m for n, m in wanted.items() if m}
    planned = evaluate(sub, wanted or None)
    sweeps = directional_sweeps(sub, max_mm=max_mm, metrics=wanted or None)
    rows = []
    for label, structs in alternatives.items():
        alt = PlanState.from_masks(state.dose, [structs[n] for n in names], state.prescription_gy)
        ev = evaluate(alt, wanted or None)
        for n in names:
            ref, new = state.masks[n], structs[n]
            c0, c1 = ref.native_cloud().centroid, new.native_cloud().centroid
            geom = {
                "volume_plan_cc": ref.mask.sum() * ref.grid.voxel_volume / 1000.0,
                "volume_alt_cc": new.mask.sum() * new.grid.voxel_volume / 1000.0,
                "centroid_shift_mm": float(np.linalg.norm(c1 - c0)),
                **{f"centroid_d{a}_mm": float(v) for a, v in zip("xyz", c1 - c0, strict=True)},
                "dice": _dice(ref, new),
            }
            geom["volume_ratio"] = geom["volume_alt_cc"] / geom["volume_plan_cc"]
            geom["rigid_comparable"] = bool(abs(geom["volume_ratio"] - 1.0) <= max_volume_change)
            for metric, value in planned.metrics[n].items():
                delta = ev.metrics[n][metric] - value
                mm, direction = equivalent_shift_mm(
                    sweeps, f"{n}:{metric}", delta, noise_floor(metric, min_delta)
                )
                if not geom["rigid_comparable"]:
                    mm, direction = float("nan"), ""  # a rigid shift cannot explain this change
                rows.append(
                    {
                        "alternative": label,
                        "structure": n,
                        **geom,
                        "metric": metric,
                        "planned": value,
                        "alternative_value": ev.metrics[n][metric],
                        "delta": delta,
                        "equivalent_shift_mm": mm,
                        "direction": direction,
                    }
                )
    return pd.DataFrame(rows)
