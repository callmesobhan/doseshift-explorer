"""Monte Carlo of fractionated treatment courses under random + systematic setup error.

One *course* = one simulated patient: a systematic error drawn once (``Sigma``) plus an
independent random error per fraction (``sigma``). Delivered dose at each structure point is the
mean over fractions of the planned dose sampled at the displaced point (van Herk's population
model under the static-dose-cloud approximation).

Two engines:

``sampled``  explicit per-fraction sampling; supports rotations and any ``n_fractions``.
``blur``     closed-form limit ``n_fractions -> inf`` for translations: planned dose convolved
             with the random-error Gaussian, then shifted by the systematic error. Fast, and a
             built-in cross-check of the sampler.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import ndimage

from .evaluate import Evaluation, evaluate, evaluate_doses
from .geometry import RigidMotion, sample_volume
from .metrics import Constraint
from .state import PlanState


@dataclass(frozen=True)
class ErrorModel:
    """Setup-error population, axes ordered (x=LR, y=AP, z=SI)."""

    sys_sigma_mm: tuple[float, float, float] = (2.0, 2.0, 2.0)  # Sigma: spread of patient means
    rand_sigma_mm: tuple[float, float, float] = (3.0, 3.0, 3.0)  # sigma: day-to-day spread
    sys_sigma_deg: tuple[float, float, float] = (0.0, 0.0, 0.0)
    rand_sigma_deg: tuple[float, float, float] = (0.0, 0.0, 0.0)
    n_fractions: int = 20

    def van_herk_margin_mm(
        self, sigma_p_mm: float = 3.2, systematic_factor: float = 2.5
    ) -> np.ndarray:
        """Margin recipe ``a Sigma + 1.64 (sqrt(sigma^2 + sp^2) - sp)`` per axis.

        ``a = 2.5`` is van Herk's 3D value (90 % of patients keep the CTV minimum dose >= 95 %);
        use ``a = 1.28`` for the 1D case in which only one axis carries error.
        """
        S, s = np.asarray(self.sys_sigma_mm), np.asarray(self.rand_sigma_mm)
        return systematic_factor * S + 1.64 * (np.sqrt(s**2 + sigma_p_mm**2) - sigma_p_mm)


@dataclass
class MCResult:
    nominal: Evaluation
    metrics: pd.DataFrame  # one row per course, columns "STRUCT:metric"
    dvh_axis: np.ndarray
    dvh_curves: dict[str, np.ndarray]  # structure -> (n_courses, len(dvh_axis))
    draws: pd.DataFrame  # systematic error per course: tx, ty, tz, rx, ry, rz
    model: ErrorModel
    method: str
    prescription_gy: float

    @property
    def n_courses(self) -> int:
        return len(self.metrics)

    def band(self, structure: str, q=(5, 50, 95)) -> np.ndarray:
        return np.percentile(self.dvh_curves[structure], q, axis=0)

    def probability(self, c: Constraint) -> float:
        """Fraction of simulated courses meeting the constraint."""
        return float(c.passes(self.metrics[c.key].to_numpy(), self.prescription_gy).mean())

    def probability_se(self, c: Constraint) -> float:
        """Binomial standard error of ``probability`` over the simulated courses."""
        p = self.probability(c)
        return float(np.sqrt(p * (1.0 - p) / self.n_courses))

    def summary(self) -> pd.DataFrame:
        nominal = self.nominal.flat()
        df = self.metrics
        out = pd.DataFrame(
            {
                "nominal": pd.Series(nominal),
                "mean": df.mean(),
                "std": df.std(),
                "p5": df.quantile(0.05),
                "p50": df.quantile(0.50),
                "p95": df.quantile(0.95),
            }
        )
        out.index.name = "structure:metric"
        return out


def _refine_axis(a: np.ndarray, axis: int, r: int) -> np.ndarray:
    """Linear node-aligned refinement of one array axis by an integer factor (exact trilinear)."""
    if r == 1:
        return a
    a = np.moveaxis(a, axis, 0)
    parts = [a[:-1] + (a[1:] - a[:-1]) * (f / r) for f in range(r)]
    out = np.concatenate([np.stack(parts, axis=1).reshape((-1,) + a.shape[1:]), a[-1:]], axis=0)
    return np.moveaxis(out, 0, axis)


def _blurred_dose(grid, data: np.ndarray, rand_sigma_mm, max_samples: float = 6e7):
    """Planned dose convolved with the random-error Gaussian, on a grid fine enough to hold it.

    * The Gaussian acts along *world* axes. The dose grid must therefore be aligned with the
      patient axes (direction cosines a signed permutation); sigma is permuted onto the grid's index
      axes. Oblique grids are rejected rather than blurred along the wrong directions.
    * If sigma is below one dose voxel on an axis, the dose is first refined (linearly, i.e. the
      very field the sampled engine interpolates) so the filter is well resolved; filtering at
      sub-voxel sigma on the original grid under-blurs badly.
    """
    sigma_world = np.asarray(rand_sigma_mm, float)
    unit = grid.affine[:3, :3] / grid.spacing_xyz  # columns: index axes (i, j, k) in world
    world_axis = np.argmax(np.abs(unit), axis=0)
    if len(set(world_axis.tolist())) != 3 or not np.allclose(
        np.abs(unit[world_axis, np.arange(3)]), 1.0, atol=1e-3
    ):
        raise ValueError(
            "the blur engine needs a dose grid aligned with the patient axes (direction cosines "
            "a signed permutation); use method='sampled' for oblique grids"
        )
    sigma_vox = sigma_world[world_axis] / grid.spacing_xyz  # per index axis (i, j, k)
    refine = np.where(
        (sigma_vox > 0) & (sigma_vox < 1.0), np.ceil(1.0 / np.maximum(sigma_vox, 1e-9)), 1
    )
    refine = refine.astype(int)  # (i, j, k)
    shape = np.array(data.shape)  # (k, j, i)
    while np.prod((shape - 1) * refine[::-1] + 1) > max_samples:
        if refine.max() == 1:
            raise ValueError("dose grid too large for the blur engine; use method='sampled'")
        refine[np.argmax(refine)] -= 1
    fine = data if refine.max() == 1 else data.astype(np.float32)  # float64 unless refined
    for array_axis in range(3):
        fine = _refine_axis(fine, array_axis, int(refine[::-1][array_axis]))
    affine = grid.affine.copy()
    affine[:3, :3] = affine[:3, :3] / refine  # divide each index-axis column by its factor
    fine_grid = type(grid)(fine.shape, affine)
    blurred = ndimage.gaussian_filter(fine, (sigma_vox * refine)[::-1], mode="constant", cval=0.0)
    return fine_grid, blurred


def simulate(
    state: PlanState,
    model: ErrorModel,
    n_courses: int = 200,
    seed: int = 0,
    method: str = "sampled",
    metrics: dict[str, list[str]] | None = None,
    dvh_axis: np.ndarray | None = None,
) -> MCResult:
    if method not in ("sampled", "blur"):
        raise ValueError("method must be 'sampled' or 'blur'")
    if method == "blur" and any(model.rand_sigma_deg):
        raise ValueError("blur engine supports random translations only")

    rng = np.random.default_rng(seed)
    names = state.names
    sizes = [len(state.clouds[n].points) for n in names]
    splits = np.cumsum(sizes)[:-1]
    all_points = np.concatenate([state.clouds[n].points for n in names])
    center = tuple(state.isocenter())
    grid, data, order = state.dose.grid, state.dose.data, state.interp_order
    top = max(1.25 * state.prescription_gy, 1.05 * float(data.max()))
    axis = np.linspace(0.0, top, 501) if dvh_axis is None else dvh_axis

    if method == "blur":
        grid, data = _blurred_dose(grid, data, model.rand_sigma_mm)

    rows, curves, draws = [], {n: [] for n in names}, []
    for _ in range(n_courses):
        t_sys = rng.normal(0.0, model.sys_sigma_mm)
        r_sys = rng.normal(0.0, model.sys_sigma_deg)
        draws.append([*t_sys, *r_sys])
        if method == "blur":
            pts = RigidMotion(tuple(t_sys), tuple(r_sys), center).apply(all_points)
            total = sample_volume(grid, data, pts, order)
        else:
            total = np.zeros(len(all_points))
            for _f in range(model.n_fractions):
                t = t_sys + rng.normal(0.0, model.rand_sigma_mm)
                r = r_sys + rng.normal(0.0, model.rand_sigma_deg)
                total += sample_volume(
                    grid, data, RigidMotion(tuple(t), tuple(r), center).apply(all_points), order
                )
            total /= model.n_fractions
        doses = dict(zip(names, np.split(total, splits), strict=True))
        ev = evaluate_doses(state, doses, metrics)
        rows.append(ev.flat())
        for n in names:
            curves[n].append(ev.dvhs[n].curve(axis))

    return MCResult(
        nominal=evaluate(state, metrics),
        metrics=pd.DataFrame(rows),
        dvh_axis=axis,
        dvh_curves={n: np.asarray(c) for n, c in curves.items()},
        draws=pd.DataFrame(draws, columns=["tx", "ty", "tz", "rx", "ry", "rz"]),
        model=model,
        method=method,
        prescription_gy=state.prescription_gy,
    )
