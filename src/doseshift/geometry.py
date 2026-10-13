"""Grid geometry and rigid patient motion.

Conventions (full write-up in docs/conventions.md):

* World space is the DICOM patient coordinate system (LPS), in mm, vectors ordered (x, y, z).
* Arrays are indexed ``[k, j, i]`` = ``[z, y, x]`` (slice, row, column), as pydicom and
  SimpleITK return them.
* Voxel ``(i, j, k)`` has its *centre* at ``affine @ [i, j, k, 1]``. DICOM
  ImagePositionPatient is the centre of the first voxel; treating it as a corner is one of the
  mistakes this project quantifies.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage
from scipy.spatial.transform import Rotation


@dataclass(frozen=True, eq=False)
class Grid:
    """A regular 3D voxel grid: array shape ``(nz, ny, nx)`` plus a 4x4 index->world affine."""

    shape: tuple[int, int, int]
    affine: np.ndarray

    def __post_init__(self) -> None:
        affine = np.asarray(self.affine, dtype=float)
        if affine.shape != (4, 4):
            raise ValueError("affine must be 4x4")
        object.__setattr__(self, "affine", affine)
        object.__setattr__(self, "shape", tuple(int(s) for s in self.shape))

    @classmethod
    def from_spacing(
        cls,
        shape: tuple[int, int, int],
        spacing_xyz: tuple[float, float, float],
        origin_xyz: tuple[float, float, float] = (0.0, 0.0, 0.0),
        direction: np.ndarray | None = None,
    ) -> Grid:
        """Axis-aligned grid unless ``direction`` (3x3, columns = unit i/j/k axes) is given."""
        d = np.eye(3) if direction is None else np.asarray(direction, dtype=float)
        affine = np.eye(4)
        affine[:3, :3] = d @ np.diag(spacing_xyz)
        affine[:3, 3] = origin_xyz
        return cls(shape, affine)

    @classmethod
    def centered(cls, shape: tuple[int, int, int], spacing_xyz: tuple[float, float, float]) -> Grid:
        """Axis-aligned grid whose voxel centres are symmetric about the world origin."""
        nz, ny, nx = shape
        origin = -0.5 * (np.array([nx, ny, nz]) - 1) * np.asarray(spacing_xyz)
        return cls.from_spacing(shape, spacing_xyz, tuple(origin))

    @property
    def spacing_xyz(self) -> np.ndarray:
        return np.linalg.norm(self.affine[:3, :3], axis=0)

    @property
    def origin_xyz(self) -> np.ndarray:
        return self.affine[:3, 3].copy()

    @property
    def voxel_volume(self) -> float:
        return float(abs(np.linalg.det(self.affine[:3, :3])))

    def index_to_world(self, ijk: np.ndarray) -> np.ndarray:
        """(N, 3) continuous index (i, j, k) -> (N, 3) world (x, y, z)."""
        ijk = np.atleast_2d(np.asarray(ijk, dtype=float))
        return ijk @ self.affine[:3, :3].T + self.affine[:3, 3]

    def world_to_index(self, xyz: np.ndarray) -> np.ndarray:
        """(N, 3) world (x, y, z) -> (N, 3) continuous index (i, j, k)."""
        xyz = np.atleast_2d(np.asarray(xyz, dtype=float))
        return (xyz - self.affine[:3, 3]) @ np.linalg.inv(self.affine[:3, :3]).T

    def mask_points(self, mask: np.ndarray) -> np.ndarray:
        """World coordinates (N, 3) of the centres of all True voxels."""
        kji = np.argwhere(mask)
        return self.index_to_world(kji[:, ::-1])

    def extent_xyz(self) -> tuple[np.ndarray, np.ndarray]:
        """Axis-aligned world bounding box of the voxel *centres*."""
        nz, ny, nx = self.shape
        corners = np.array(
            [[i, j, k] for i in (0, nx - 1) for j in (0, ny - 1) for k in (0, nz - 1)], float
        )
        w = self.index_to_world(corners)
        return w.min(axis=0), w.max(axis=0)

    def translated_voxels(self, voxels_ijk: tuple[float, float, float]) -> Grid:
        """Same array, but located ``voxels_ijk`` (index units) further along the grid axes."""
        affine = self.affine.copy()
        affine[:3, 3] += affine[:3, :3] @ np.asarray(voxels_ijk, dtype=float)
        return Grid(self.shape, affine)

    def allclose(self, other: Grid, atol: float = 1e-6) -> bool:
        return self.shape == other.shape and np.allclose(self.affine, other.affine, atol=atol)


@dataclass(frozen=True)
class RigidMotion:
    """Rigid *patient* motion: anatomy planned at ``p`` is found at ``R (p - c) + c + t``.

    ``rotation_deg`` are right-handed rotations about the fixed world x, y, z axes, applied in
    that order (``R = Rz Ry Rx``). A positive ``translation_mm[0]`` moves the patient towards
    the patient's left (+x in LPS). The dose distribution stays fixed in the room.
    """

    translation_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)
    rotation_deg: tuple[float, float, float] = (0.0, 0.0, 0.0)
    center_mm: tuple[float, float, float] = (0.0, 0.0, 0.0)

    @property
    def matrix(self) -> np.ndarray:
        return Rotation.from_euler("xyz", self.rotation_deg, degrees=True).as_matrix()

    @classmethod
    def from_matrix(cls, R: np.ndarray, translation_mm, center_mm=(0.0, 0.0, 0.0)) -> RigidMotion:
        euler = Rotation.from_matrix(R).as_euler("xyz", degrees=True)
        return cls(tuple(translation_mm), tuple(euler), tuple(center_mm))

    @property
    def is_identity(self) -> bool:
        return not (any(self.translation_mm) or any(self.rotation_deg))

    def apply(self, points: np.ndarray) -> np.ndarray:
        c = np.asarray(self.center_mm, dtype=float)
        t = np.asarray(self.translation_mm, dtype=float)
        return (points - c) @ self.matrix.T + c + t

    def inverse_apply(self, points: np.ndarray) -> np.ndarray:
        c = np.asarray(self.center_mm, dtype=float)
        t = np.asarray(self.translation_mm, dtype=float)
        return (points - c - t) @ self.matrix + c

    def warp_volume(
        self, grid: Grid, array: np.ndarray, order: int = 1, cval: float = 0.0
    ) -> np.ndarray:
        """Image of ``array`` after the patient moved: ``out(q) = array(motion^-1(q))``."""
        nz, ny, nx = grid.shape
        kk, jj, ii = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
        q = grid.index_to_world(np.stack([ii.ravel(), jj.ravel(), kk.ravel()], axis=1))
        ijk = grid.world_to_index(self.inverse_apply(q))
        out = ndimage.map_coordinates(
            array, [ijk[:, 2], ijk[:, 1], ijk[:, 0]], order=order, mode="constant", cval=cval
        )
        return out.reshape(grid.shape)


def sample_volume(
    grid: Grid, array: np.ndarray, points_xyz: np.ndarray, order: int = 1, cval: float = 0.0
) -> np.ndarray:
    """Interpolate ``array`` at world points; ``cval`` (0 Gy for dose) outside the grid."""
    ijk = grid.world_to_index(points_xyz)
    return ndimage.map_coordinates(
        array, [ijk[:, 2], ijk[:, 1], ijk[:, 0]], order=order, mode="constant", cval=cval
    )
