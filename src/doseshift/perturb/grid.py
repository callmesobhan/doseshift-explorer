"""Numerical perturbations: alternative readings of the dose/structure grids.

None of these move the patient. They change how dose values are *located* or *sampled*, and are
the "dose-grid alignment assumptions" this project quantifies.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from ..geometry import Grid
from ..state import DoseGrid, MaskStructure, PlanState, PointCloud
from .base import Perturbation


@dataclass(frozen=True)
class DoseGridOffset(Perturbation):
    """Locate the dose array ``voxels`` (index units along i/j/k) further along its own axes.

    ``+0.5`` on each axis is what you get by treating DICOM ImagePositionPatient (a voxel
    *centre*) as a voxel corner.

    In the room this displaces the dose by ``A @ v``, where ``A`` is the 3x3 part of the grid
    affine (so it follows the grid's orientation and slice direction). It is equivalent to
    moving the patient by ``t = -A @ v``. Only for an axis-aligned grid with positive steps
    does that reduce to ``-v * spacing`` along x, y, z.
    """

    voxels: tuple[float, float, float] = (0.5, 0.5, 0.5)
    kind = "numerical"

    def apply(self, state: PlanState) -> PlanState:
        dose = DoseGrid(state.dose.grid.translated_voxels(self.voxels), state.dose.data)
        return state.evolve(self.describe(), dose=dose)

    def describe(self) -> str:
        return f"dose-offset[{','.join(f'{v:+g}' for v in self.voxels)}]vox(ijk)"


@dataclass(frozen=True)
class DoseCoarsen(Perturbation):
    """Block-average the dose onto a grid ``factor`` times coarser (lower calc resolution)."""

    factor: int = 2
    kind = "numerical"

    def apply(self, state: PlanState) -> PlanState:
        f = int(self.factor)
        grid, data = state.dose.grid, state.dose.data
        nz, ny, nx = (s // f * f for s in data.shape)
        if min(nz, ny, nx) == 0:
            raise ValueError("coarsening factor larger than the dose grid")
        if (nz, ny, nx) != data.shape:
            warnings.warn(
                f"dose shape {data.shape} is not a multiple of {f}; trailing voxels are dropped, "
                "so structures near that edge read 0 Gy (an artefact of coarsening)",
                stacklevel=1,
            )
        blocks = data[:nz, :ny, :nx].reshape(nz // f, f, ny // f, f, nx // f, f)
        coarse = blocks.mean(axis=(1, 3, 5))
        affine = grid.affine.copy()
        affine[:3, 3] += affine[:3, :3] @ np.full(3, (f - 1) / 2.0)  # centre of first block
        affine[:3, :3] = affine[:3, :3] * f
        new = DoseGrid(Grid(coarse.shape, affine), coarse)
        return state.evolve(self.describe(), dose=new)

    def describe(self) -> str:
        return f"dose-coarsen x{self.factor}"


@dataclass(frozen=True)
class InterpolationOrder(Perturbation):
    """Dose interpolation when sampling at structure points: 0 nearest, 1 trilinear, 3 cubic."""

    order: int = 1
    kind = "numerical"

    def apply(self, state: PlanState) -> PlanState:
        return state.evolve(self.describe(), interp_order=int(self.order))

    def describe(self) -> str:
        return f"interp order={self.order}"


@dataclass(frozen=True)
class StructureSampling(Perturbation):
    """Evaluate structures on the *dose* grid instead of the structure (CT) grid.

    * ``dose_center``     - binary: a dose voxel is in the structure if its centre is.
    * ``dose_fractional`` - each dose voxel is weighted by the fraction of it inside the
      structure, estimated from ``subvoxel``^3 sub-samples.

    Must precede any motion perturbation: it rebuilds the clouds from the reference masks.
    """

    mode: str = "dose_center"
    subvoxel: int = 3
    kind = "numerical"

    def apply(self, state: PlanState) -> PlanState:
        if state.motion_applied:
            raise ValueError("StructureSampling must be applied before any motion perturbation")
        if self.mode not in ("dose_center", "dose_fractional"):
            raise ValueError(f"unknown mode {self.mode!r}")
        s = self.subvoxel if self.mode == "dose_fractional" else 1
        clouds = {n: _cloud_on_grid(m, state.dose.grid, s) for n, m in state.masks.items()}
        return state.evolve(self.describe(), clouds=clouds)

    def describe(self) -> str:
        return f"structures-on-dose-grid [{self.mode}]"


def mask_contains(struct: MaskStructure, points: np.ndarray) -> np.ndarray:
    idx = np.rint(struct.grid.world_to_index(points)).astype(int)
    shape_ijk = np.array(struct.grid.shape[::-1])
    inside = np.all((idx >= 0) & (idx < shape_ijk), axis=1)
    out = np.zeros(len(points), dtype=bool)
    i, j, k = idx[inside].T
    out[inside] = struct.mask[k, j, i]
    return out


def _cloud_on_grid(struct: MaskStructure, grid: Grid, sub: int) -> PointCloud:
    """Structure as weighted points at dose-voxel centres.

    Candidate voxels are *not* clipped to the dose grid: the part of a structure outside it stays
    in the cloud (and reads 0 Gy exactly as the native evaluation does), so the structure's volume
    is conserved and the perturbation changes sampling only, never the structure definition.
    """
    pts = struct.grid.mask_points(struct.mask)
    ijk = grid.world_to_index(pts)
    # a structure voxel's reach along each dose index axis sets the margin around the point cloud
    step = np.linalg.inv(grid.affine[:3, :3]) @ struct.grid.affine[:3, :3]
    margin = np.ceil(0.5 * np.abs(step).sum(axis=1)).astype(int) + 1
    lo = np.floor(ijk.min(axis=0)).astype(int) - margin
    hi = np.ceil(ijk.max(axis=0)).astype(int) + margin
    axes = [np.arange(a, b + 1) for a, b in zip(lo, hi, strict=True)]
    cand = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3).astype(float)
    off = (np.arange(sub) + 0.5) / sub - 0.5
    offsets = np.stack(np.meshgrid(off, off, off, indexing="ij"), axis=-1).reshape(-1, 3)
    occupancy = np.zeros(len(cand))
    for o in offsets:
        occupancy += mask_contains(struct, grid.index_to_world(cand + o))
    occupancy /= len(offsets)
    keep = occupancy > 0
    return PointCloud(grid.index_to_world(cand[keep]), occupancy[keep] * grid.voxel_volume)
