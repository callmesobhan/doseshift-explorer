"""Physical perturbations: the patient is not where the plan assumed."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage

from ..geometry import Grid, RigidMotion
from ..state import PlanState, PointCloud
from .base import Perturbation


def _moved(state: PlanState, fn, note: str) -> PlanState:
    clouds = {n: PointCloud(fn(c.points), c.volumes) for n, c in state.clouds.items()}
    return state.evolve(note, clouds=clouds, motion_applied=True)


@dataclass(frozen=True)
class RigidShift(Perturbation):
    """Translate and/or rotate the patient about ``center_mm`` (default: PTV centroid)."""

    translation_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)
    rotation_deg: tuple[float, float, float] = (0.0, 0.0, 0.0)
    center_mm: tuple[float, float, float] | None = None
    kind = "physical"

    def motion_for(self, state: PlanState) -> RigidMotion:
        center = self.center_mm if self.center_mm is not None else tuple(state.isocenter())
        return RigidMotion(self.translation_mm, self.rotation_deg, center)

    def apply(self, state: PlanState) -> PlanState:
        motion = self.motion_for(state)
        return _moved(state, motion.apply, self.describe())

    def describe(self) -> str:
        t = ",".join(f"{v:g}" for v in self.translation_mm)
        r = ",".join(f"{v:g}" for v in self.rotation_deg)
        return f"shift[{t}]mm rot[{r}]deg"


@dataclass(frozen=True)
class SmoothDeformation(Perturbation):
    """Random smooth, zero-mean displacement field; realised RMS ``rms_mm`` per component.

    Correlation length ``corr_mm``. The common translation is removed (rigid motion is a separate
    perturbation) and the RMS is measured over all structure points. A controlled stand-in for
    anatomical change (filling, weight loss) when no deformable registration is available.
    Deterministic for a given ``seed``.
    """

    rms_mm: float = 2.0
    corr_mm: float = 25.0
    seed: int = 0
    kind = "physical"

    def apply(self, state: PlanState) -> PlanState:
        allp = np.concatenate([c.points for c in state.clouds.values()])
        h = max(self.corr_mm / 3.0, 1.0)
        lo, hi = allp.min(axis=0) - 2 * self.corr_mm, allp.max(axis=0) + 2 * self.corr_mm
        n = np.maximum(np.ceil((hi - lo) / h).astype(int) + 1, 4)
        field_grid = Grid.from_spacing((int(n[2]), int(n[1]), int(n[0])), (h, h, h), tuple(lo))
        rng = np.random.default_rng(self.seed)
        comps = [
            ndimage.gaussian_filter(rng.standard_normal(field_grid.shape), self.corr_mm / h)
            for _ in range(3)
        ]

        def raw(points: np.ndarray) -> np.ndarray:
            ijk = field_grid.world_to_index(points)
            coords = [ijk[:, 2], ijk[:, 1], ijk[:, 0]]
            return np.stack(
                [ndimage.map_coordinates(c, coords, order=1, mode="nearest") for c in comps],
                axis=1,
            )

        # Normalise over the points that are actually displaced: remove the common (rigid)
        # translation so this is purely non-rigid, then scale so the realised RMS displacement
        # per component equals ``rms_mm``. (Normalising the whole field instead gives a realised
        # RMS far from the request when ``corr_mm`` is large relative to the structures.)
        d_all = raw(allp)
        mean = d_all.mean(axis=0)
        scale = self.rms_mm / max(float(np.sqrt(((d_all - mean) ** 2).mean())), 1e-12)

        def displace(points: np.ndarray) -> np.ndarray:
            return points + (raw(points) - mean) * scale

        return _moved(state, displace, self.describe())

    def describe(self) -> str:
        return f"deform[rms={self.rms_mm:g}mm corr={self.corr_mm:g}mm seed={self.seed}]"
