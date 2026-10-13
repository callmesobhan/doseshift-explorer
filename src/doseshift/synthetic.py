"""Analytic prostate-like phantom: runs everywhere, needs no patient data.

The dose is a smooth plateau with an erfc penumbra around the PTV, so every DVH feature has a
known cause. CT and structure masks live on a 2x2x3 mm grid; the dose lives on a *different*
3 mm isotropic grid (voxel centres do not coincide), as in clinical exports.

LPS axes: +x patient left, +y posterior, +z superior.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage
from scipy.special import ndtr, ndtri

from .geometry import Grid
from .state import DoseGrid, MaskStructure, PlanState

PRESCRIPTION_GY = 60.0
PENUMBRA_SIGMA_MM = 3.2  # van Herk's photon penumbra sigma
PLATEAU = 1.04  # plateau dose / Rx
COVERAGE_AT_SURFACE = 0.95  # 95 % isodose on the PTV surface


@dataclass(frozen=True)
class Ellipsoid:
    center: tuple[float, float, float]
    radii: tuple[float, float, float]

    def rho(self, xyz: np.ndarray) -> np.ndarray:
        return np.sqrt((((xyz - np.asarray(self.center)) / np.asarray(self.radii)) ** 2).sum(-1))

    def contains(self, xyz: np.ndarray) -> np.ndarray:
        return self.rho(xyz) <= 1.0

    def signed_distance(self, xyz: np.ndarray) -> np.ndarray:
        """First-order distance to the surface in mm (negative inside)."""
        r, c = np.asarray(self.radii), np.asarray(self.center)
        rho = np.maximum(self.rho(xyz), 1e-6)
        grad = np.linalg.norm((xyz - c) / r**2, axis=-1) / rho
        sd = (rho - 1.0) / np.maximum(grad, 1e-9)
        return np.maximum(sd, -min(self.radii))  # finite at the exact centre

    def shrunk(self, mm: float) -> Ellipsoid:
        return Ellipsoid(self.center, tuple(r - mm for r in self.radii))


PTV = Ellipsoid((0.0, 0.0, 0.0), (28.0, 24.0, 26.0))
RECTUM = Ellipsoid((0.0, 41.0, 0.0), (18.0, 14.0, 40.0))  # 3 mm behind the PTV surface
BLADDER = Ellipsoid((0.0, -49.0, 0.0), (26.0, 20.0, 24.0))  # 5 mm in front


@dataclass
class PhantomCase:
    ct_grid: Grid
    ct_hu: np.ndarray
    dose: DoseGrid
    masks: dict[str, MaskStructure]
    prescription_gy: float

    def to_state(self) -> PlanState:
        return PlanState.from_masks(self.dose, list(self.masks.values()), self.prescription_gy)


def _coords(grid: Grid) -> np.ndarray:
    nz, ny, nx = grid.shape
    kk, jj, ii = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    ijk = np.stack([ii, jj, kk], axis=-1).reshape(-1, 3)
    return grid.index_to_world(ijk).reshape(nz, ny, nx, 3)


def make_dose(grid: Grid, seed: int = 0, noise: float = 0.01) -> np.ndarray:
    """PLATEAU * Rx inside, erfc penumbra outside, 95 % isodose on the PTV surface."""
    s0 = ndtri(COVERAGE_AT_SURFACE / PLATEAU) * PENUMBRA_SIGMA_MM
    sd = PTV.signed_distance(_coords(grid))
    dose = PRESCRIPTION_GY * PLATEAU * ndtr((s0 - sd) / PENUMBRA_SIGMA_MM)
    if noise:
        rng = np.random.default_rng(seed)
        f = ndimage.gaussian_filter(rng.standard_normal(grid.shape), 2.0)
        dose = dose * (1.0 + noise * f / f.std())
    return dose


def make_ct(grid: Grid, seed: int = 0) -> np.ndarray:
    """Pelvis-like CT with enough texture for rigid registration to lock onto."""
    xyz = _coords(grid)
    body = Ellipsoid((0, 0, 0), (80, 78, 66))
    bone = [Ellipsoid((s * 58.0, 5.0, 0.0), (20, 20, 20)) for s in (-1, 1)]
    bone.append(Ellipsoid((0.0, 72.0, 0.0), (24, 7, 45)))
    rng = np.random.default_rng(seed)
    tex = ndimage.gaussian_filter(rng.standard_normal(grid.shape), 3.0)
    hu = np.full(grid.shape, -1000.0)
    hu[body.contains(xyz)] = 30.0
    hu[PTV.contains(xyz)] = 45.0
    for b in bone:
        hu[b.contains(xyz)] = 650.0
    hu[body.contains(xyz)] += 60.0 * tex[body.contains(xyz)] / tex.std()
    return ndimage.gaussian_filter(hu, 0.8)


def make_case(seed: int = 0, ctv_margin_mm: float = 5.0) -> PhantomCase:
    ct_grid = Grid.centered((48, 90, 90), (2.0, 2.0, 3.0))
    dose_grid = Grid.centered((48, 60, 60), (3.0, 3.0, 3.0))
    xyz = _coords(ct_grid)
    ptv = PTV.contains(xyz)
    masks = {
        "PTV": MaskStructure("PTV", "ptv", ptv, ct_grid),
        "CTV": MaskStructure("CTV", "ctv", PTV.shrunk(ctv_margin_mm).contains(xyz), ct_grid),
        "Rectum": MaskStructure("Rectum", "oar", RECTUM.contains(xyz) & ~ptv, ct_grid),
        "Bladder": MaskStructure("Bladder", "oar", BLADDER.contains(xyz) & ~ptv, ct_grid),
    }
    return PhantomCase(
        ct_grid,
        make_ct(ct_grid, seed),
        DoseGrid(dose_grid, make_dose(dose_grid, seed)),
        masks,
        PRESCRIPTION_GY,
    )


def make_state(seed: int = 0, ctv_margin_mm: float = 5.0) -> PlanState:
    return make_case(seed, ctv_margin_mm).to_state()
