"""Immutable plan state that perturbations transform and evaluators read.

``PlanState`` separates three things so that physical and numerical perturbations stay
orthogonal:

* ``dose``   - the planned dose on its own grid (room frame, never moved by patient motion);
* ``masks``  - reference structures on their own grid (kept so sampling can be rebuilt);
* ``clouds`` - structures as weighted point clouds, the thing patient motion acts on. Dose is
  evaluated by interpolating the dose grid at the cloud points (static-dose-cloud approximation).
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field, replace

import numpy as np

from .geometry import Grid, sample_volume

TARGET_ROLES = ("ptv", "ctv")
ROLES = (*TARGET_ROLES, "oar")


@dataclass(frozen=True, eq=False)
class DoseGrid:
    grid: Grid
    data: np.ndarray  # Gy, shape == grid.shape

    def __post_init__(self) -> None:
        if tuple(self.data.shape) != self.grid.shape:
            raise ValueError(f"dose shape {self.data.shape} != grid shape {self.grid.shape}")


@dataclass(frozen=True, eq=False)
class MaskStructure:
    name: str
    role: str
    mask: np.ndarray  # bool, shape == grid.shape
    grid: Grid

    def __post_init__(self) -> None:
        if self.role not in ROLES:
            raise ValueError(f"role must be one of {ROLES}, got {self.role!r}")
        if tuple(self.mask.shape) != self.grid.shape:
            raise ValueError(f"mask shape {self.mask.shape} != grid shape {self.grid.shape}")

    def native_cloud(self) -> PointCloud:
        return PointCloud(
            self.grid.mask_points(self.mask), np.full(int(self.mask.sum()), self.grid.voxel_volume)
        )


@dataclass(frozen=True, eq=False)
class PointCloud:
    points: np.ndarray  # (N, 3) world mm
    volumes: np.ndarray  # (N,) mm^3 represented by each point

    @property
    def total_volume(self) -> float:
        return float(self.volumes.sum())

    @property
    def centroid(self) -> np.ndarray:
        return (self.points * self.volumes[:, None]).sum(axis=0) / self.volumes.sum()


@dataclass(frozen=True, eq=False)
class PlanState:
    dose: DoseGrid
    masks: dict[str, MaskStructure]
    clouds: dict[str, PointCloud]
    prescription_gy: float
    rx_label: str = "Rx"  # "reference dose" when the data holds no prescription
    interp_order: int = 1
    motion_applied: bool = False
    history: tuple[str, ...] = field(default_factory=tuple)

    @classmethod
    def from_masks(
        cls, dose: DoseGrid, structures: list[MaskStructure], prescription_gy: float
    ) -> PlanState:
        if not (math.isfinite(prescription_gy) and prescription_gy > 0):
            raise ValueError(f"prescription_gy must be finite and > 0, got {prescription_gy}")
        masks = {s.name: s for s in structures}
        if len(masks) != len(structures):
            raise ValueError("structure names must be unique")
        clouds = {}
        for n, s in masks.items():
            if not s.mask.any():
                raise ValueError(f"structure {n!r} is empty (no voxels inside its contours)")
            clouds[n] = s.native_cloud()
        state = cls(dose, masks, clouds, float(prescription_gy))
        for n in masks:
            frac = state.outside_dose_fraction(n)
            if frac > 0.005:
                warnings.warn(
                    f"{100 * frac:.1f} % of structure {n!r} lies outside the dose grid and is "
                    "evaluated as 0 Gy; its DVH will be wrong",
                    stacklevel=2,
                )
        return state

    def outside_dose_fraction(self, name: str) -> float:
        """Volume fraction of a structure outside the dose grid (those points read 0 Gy)."""
        c = self.clouds[name]
        ijk = self.dose.grid.world_to_index(c.points)
        upper = np.array(self.dose.grid.shape[::-1]) - 1
        inside = np.all((ijk >= 0) & (ijk <= upper), axis=1)
        return float(c.volumes[~inside].sum() / c.volumes.sum())

    def evolve(self, note: str | None = None, **changes) -> PlanState:
        history = self.history + ((note,) if note else ())
        return replace(self, history=history, **changes)

    @property
    def names(self) -> list[str]:
        return list(self.masks)

    def role(self, name: str) -> str:
        return self.masks[name].role

    def isocenter(self) -> np.ndarray:
        """Centroid of the PTV (else any target) in the planned position."""
        for roles in (("ptv",), TARGET_ROLES):
            pts = [self.clouds[n] for n, m in self.masks.items() if m.role in roles]
            if pts:
                allp = np.concatenate([c.points for c in pts])
                allv = np.concatenate([c.volumes for c in pts])
                return (allp * allv[:, None]).sum(axis=0) / allv.sum()
        return np.asarray(self.dose.grid.extent_xyz()).mean(axis=0)

    def dose_at(self, points: np.ndarray) -> np.ndarray:
        return sample_volume(self.dose.grid, self.dose.data, points, self.interp_order)

    def structure_doses(self) -> dict[str, np.ndarray]:
        return {n: self.dose_at(c.points) for n, c in self.clouds.items()}
