"""Turn a (perturbed) PlanState, or externally computed doses, into DVHs and metrics."""

from __future__ import annotations

import warnings
import weakref
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .dvh import DVH
from .metrics import DEFAULT_METRICS, MetricContext, compute_metric
from .state import PlanState


@dataclass
class Evaluation:
    dvhs: dict[str, DVH]
    metrics: dict[str, dict[str, float]]  # structure -> metric -> value
    history: tuple[str, ...] = field(default_factory=tuple)

    def flat(self) -> dict[str, float]:
        """``{"PTV:D95%": 58.1, ...}`` - the row format used by sweeps and Monte Carlo."""
        return {f"{s}:{m}": v for s, ms in self.metrics.items() for m, v in ms.items()}

    def to_frame(self) -> pd.DataFrame:
        rows = [(s, m, v) for s, ms in self.metrics.items() for m, v in ms.items()]
        return pd.DataFrame(rows, columns=["structure", "metric", "value"])


MIN_SAMPLES = 30


def resolve_metrics(state: PlanState, metrics: dict[str, list[str]] | None) -> dict[str, list[str]]:
    """Per-structure metric lists: explicit entry by name, else the default for its role."""
    metrics = metrics or {}
    out = {}
    for name in state.names:
        out[name] = list(metrics.get(name, DEFAULT_METRICS[state.role(name)]))
    return out


_ISODOSE_CACHE: weakref.WeakKeyDictionary = weakref.WeakKeyDictionary()


def _refine_axis(a: np.ndarray, axis: int, r: int) -> np.ndarray:
    """Linear (node-aligned) refinement of one array axis by an integer factor."""
    if r == 1:
        return a
    a = np.moveaxis(a, axis, 0)
    lo, hi = a[:-1], a[1:]
    parts = [lo + (hi - lo) * (f / r) for f in range(r)]
    out = np.stack(parts, axis=1).reshape((-1,) + a.shape[1:])
    out = np.concatenate([out, a[-1:]], axis=0)
    return np.moveaxis(out, 0, axis)


def _isodose_volumes(
    state: PlanState, target_step_mm: float = 1.0, max_samples: float = 6e7
) -> tuple[float, float]:
    grid, data = state.dose.grid, state.dose.data
    rx = state.prescription_gy
    inside = data >= 0.5 * rx
    if not inside.any():
        return 0.0, 0.0
    idx = np.argwhere(inside)  # array axes (k, j, i)
    lo = np.maximum(idx.min(axis=0) - 1, 0)
    hi = np.minimum(idx.max(axis=0) + 1, np.array(data.shape) - 1)
    box = data[lo[0] : hi[0] + 1, lo[1] : hi[1] + 1, lo[2] : hi[2] + 1]
    spacing = grid.spacing_xyz[::-1]  # (k, j, i)
    r = np.maximum(1, np.ceil(spacing / target_step_mm)).astype(int)
    if state.interp_order == 0:
        r = np.ones(3, dtype=int)  # nearest: the nodes already are the field
    while np.prod((np.array(box.shape) - 1) * r + 1) > max_samples and r.max() > 1:
        r[np.argmax(r)] -= 1
    fine = box.astype(np.float32)
    for axis in range(3):
        fine = _refine_axis(fine, axis, int(r[axis]))
    cell = float(np.prod(spacing / r))
    return float((fine >= rx).sum() * cell), float((fine >= 0.5 * rx).sum() * cell)


def isodose_volumes_mm3(state: PlanState) -> tuple[float, float]:
    """(PIV, V50): volume of the *interpolated* dose field >= Rx and >= 0.5 Rx.

    Evaluated on a lattice refined to about 1 mm so that it is the same field (and the same
    interpolation) that ``TV_PIV`` samples at structure points; counting dose-grid nodes instead
    biases Paddick's CI with dose-grid spacing. Independent of patient motion; cached per dose grid.
    """
    cache = _ISODOSE_CACHE.setdefault(state.dose, {})
    key = (state.prescription_gy, state.interp_order)
    if key not in cache:
        cache[key] = _isodose_volumes(state)
    return cache[key]


def evaluate_doses(
    state: PlanState, doses: dict[str, np.ndarray], metrics: dict[str, list[str]] | None = None
) -> Evaluation:
    """Evaluate structure doses sampled by the caller (used directly by the Monte Carlo loop)."""
    wanted = resolve_metrics(state, metrics)
    piv, v50 = isodose_volumes_mm3(state)
    dvhs, results = {}, {}
    for name, d in doses.items():
        vols = state.clouds[name].volumes
        if len(vols) == 0:
            raise ValueError(
                f"structure {name!r} has no sample points (too small for this sampling?)"
            )
        if len(vols) < MIN_SAMPLES:
            warnings.warn(
                f"structure {name!r} has only {len(vols)} sample points; its DVH metrics are "
                "dominated by sampling (a change here may measure quantisation, not dose)",
                stacklevel=2,
            )
        dvh = DVH(d, vols)
        ctx = MetricContext(dvh, d, vols, state.prescription_gy, piv, v50)
        dvhs[name] = dvh
        results[name] = {m: compute_metric(m, ctx) for m in wanted[name]}
    return Evaluation(dvhs, results, state.history)


def evaluate(state: PlanState, metrics: dict[str, list[str]] | None = None) -> Evaluation:
    return evaluate_doses(state, state.structure_doses(), metrics)
