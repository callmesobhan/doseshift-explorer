"""Shared test helpers (import as ``from helpers import ...``)."""

import numpy as np


def cbct_like(fixed_grid, values, motion, offset, seed):
    """The same anatomy as another scanner would show it: coarser anisotropic grid, an unrelated
    frame origin (``offset``), a nonlinear intensity remap, blur and noise."""
    from scipy import ndimage

    from doseshift.geometry import Grid, sample_volume

    moved = motion.warp_volume(fixed_grid, values, order=1, cval=-1000.0)  # patient moved
    lo, hi = fixed_grid.extent_xyz()
    sp = (3.0, 3.0, 4.0)
    n = tuple(np.ceil((hi - lo) / sp).astype(int) + 1)
    g = Grid.from_spacing((n[2], n[1], n[0]), sp, tuple(lo + offset))
    kk, jj, ii = np.meshgrid(np.arange(n[2]), np.arange(n[1]), np.arange(n[0]), indexing="ij")
    q = g.index_to_world(np.stack([ii.ravel(), jj.ravel(), kk.ravel()], 1))
    v = sample_volume(fixed_grid, moved, q - offset, order=1, cval=-1000.0).reshape(g.shape)
    v = 1500.0 * (np.clip(v + 1000.0, 0, None) / 4000.0) ** 0.55 - 1000.0
    rng = np.random.default_rng(seed)
    return g, (ndimage.gaussian_filter(v, 0.7) + rng.normal(0, 25, v.shape)).astype(np.float32)
