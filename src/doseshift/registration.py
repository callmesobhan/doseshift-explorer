"""Rigid image registration (SimpleITK) expressed in the project's motion convention.

Registering a *treatment* image (moving) to the *planning* CT (fixed) yields a transform that
maps planned positions to where the same anatomy sits on the day: exactly ``RigidMotion``. That
lets measured setup errors feed the Monte Carlo error model directly, and lets the perturbation
engine be validated against an independent implementation (tests/test_registration.py).

Requires the ``registration`` extra: ``pip install doseshift[registration]``.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

from .geometry import Grid, RigidMotion

try:
    import SimpleITK as sitk
except ImportError:  # pragma: no cover - exercised only without the extra
    sitk = None


def _require_sitk():
    if sitk is None:
        raise ImportError("SimpleITK is required: pip install 'doseshift[registration]'")
    return sitk


def to_sitk(grid: Grid, array: np.ndarray):
    """Array [k, j, i] on ``grid`` -> sitk.Image with matching origin / spacing / direction."""
    s = _require_sitk()
    img = s.GetImageFromArray(np.ascontiguousarray(array, dtype=np.float32))
    spacing = grid.spacing_xyz
    img.SetSpacing(tuple(float(v) for v in spacing))
    img.SetOrigin(tuple(float(v) for v in grid.origin_xyz))
    img.SetDirection(tuple(float(v) for v in (grid.affine[:3, :3] / spacing).ravel()))
    return img


@dataclass(frozen=True)
class RegistrationResult:
    """Outcome of ``register_rigid``.

    ``motion`` is only meaningful if ``converged``; ``iterations`` / ``stop_condition`` describe the
    finest level, ``levels`` lists (iterations, stop condition) for every resolution level, and
    ``overlap_fraction`` is the share of the fixed image's foreground that falls inside the moving
    image at the starting transform.
    """

    motion: RigidMotion
    metric_value: float
    iterations: int
    stop_condition: str
    converged: bool = True
    overlap_fraction: float = 1.0
    levels: tuple[tuple[int, str], ...] = ()


def otsu_threshold(values: np.ndarray, bins: int = 256) -> float:
    """Otsu's threshold: the intensity that best separates a dark background from the body."""
    hist, edges = np.histogram(values[np.isfinite(values)], bins=bins)
    centres = 0.5 * (edges[:-1] + edges[1:])
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * centres) / np.maximum(w0, 1)
    m1 = (np.sum(hist * centres) - np.cumsum(hist * centres)) / np.maximum(w1, 1)
    between = w0 * w1 * (m0 - m1) ** 2
    return float(centres[int(np.argmax(between))])


def intensity_centroid_mm(
    grid: Grid, values: np.ndarray, threshold: float | None = None
) -> np.ndarray:
    """World position of the centre of mass of the voxels above ``threshold``.

    ``threshold=None`` uses Otsu's method, so it works for any intensity scale (HU, a CBCT's
    different scaling, MR). Use it as an initial translation between images whose coordinate
    frames are unrelated (e.g. a CBCT centred on the isocentre and a planning CT with its own
    origin). It is only a starting point: a truncated field of view biases it.
    """
    if threshold is None:
        threshold = otsu_threshold(values)
    return grid.index_to_world(np.argwhere(values > threshold).mean(axis=0)[::-1][None])[0]


def _overlap_fraction(fixed, moving, translation) -> float:
    """Share of the fixed image's foreground that lands inside the moving image's extent when the
    fixed points are displaced by ``translation`` (rotation ignored: it is a start-up check)."""
    fg, fv = fixed
    pts = fg.index_to_world(np.argwhere(fv > otsu_threshold(fv))[::7, ::-1]) + np.asarray(
        translation
    )
    mg = moving[0]
    ijk = mg.world_to_index(pts)
    upper = np.array(mg.shape[::-1]) - 1
    return float(np.mean(np.all((ijk >= 0) & (ijk <= upper), axis=1))) if len(pts) else 0.0


def _usable_levels(shape, shrink_factors, sigmas):
    """Drop pyramid levels that would shrink the smallest image dimension below ~4 voxels."""
    limit = max(1, min(shape) // 4)
    keep = [(f, g) for f, g in zip(shrink_factors, sigmas, strict=True) if f <= limit]
    return keep or [(1, 0.0)]


def register_rigid(
    fixed: tuple[Grid, np.ndarray],
    moving: tuple[Grid, np.ndarray],
    center_mm: tuple[float, float, float] | None = None,
    with_rotation: bool = True,
    shrink_factors: tuple[int, ...] = (4, 2, 1),
    smoothing_sigmas_vox: tuple[float, ...] = (2.0, 1.0, 0.0),
    metric: str = "mse",
    initial_translation_mm: tuple[float, float, float] | None = None,
    sampling_percentage: float = 0.25,
    seed: int = 0,
) -> RegistrationResult:
    """Rigid registration of ``moving`` onto ``fixed``; returns the motion of the anatomy.

    The result maps a point of the fixed image to the same anatomy in the moving image,
    ``p_moving = R (p - c) + c + t``. Fixed and moving may sit on different grids and in different
    coordinate frames; then ``t`` contains the frame offset and is *not* a setup error unless the
    frames are related (give ``initial_translation_mm`` to start near the answer).

    ``metric="mse"``  mean squares: same modality and intensity scale (CT vs CT, MR vs MR).
    ``metric="mi"``   Mattes mutual information with seeded random sampling: different scanners or
                      intensity scales (CT vs CBCT). Slower and noisier; check ``stop_condition``.

    ``center_mm`` is the rotation centre in the *fixed* frame (use the PTV centroid so rotation and
    translation are reported about the isocentre); default is the fixed image centre. The reported
    translation is the displacement *of that centre*, so it is only a "shift at the isocentre"
    if you pass the isocentre; the motion itself (which carries its centre) is correct either way.

    Raises ``ValueError`` if the images barely overlap at the start (unrelated frames) and
    ``RuntimeError`` if the metric was never evaluated; warns if a level hit the iteration limit
    (``converged`` is then False).
    """
    if metric not in ("mse", "mi"):
        raise ValueError("metric must be 'mse' or 'mi'")
    for name, (grid, vals) in (("fixed", fixed), ("moving", moving)):
        if not np.isfinite(vals).all():
            raise ValueError(f"{name} image contains NaN or infinite values")
        if min(grid.shape) < 4:
            raise ValueError(
                f"{name} image has only {min(grid.shape)} voxels along one axis; ITK needs at "
                "least 4 along every axis to filter and register a volume"
            )
    s = _require_sitk()
    f_img, m_img = to_sitk(*fixed), to_sitk(*moving)
    if center_mm is None:
        lo, hi = fixed[0].extent_xyz()
        center_mm = tuple((lo + hi) / 2)
    start = np.zeros(3) if initial_translation_mm is None else np.asarray(initial_translation_mm)
    overlap = _overlap_fraction(fixed, moving, start)
    if overlap < 0.2:
        raise ValueError(
            f"fixed and moving images barely overlap at the starting transform ({overlap:.0%} of "
            "the fixed foreground lies inside the moving image): their coordinate frames are "
            "probably unrelated. Pass initial_translation_mm, e.g. "
            "intensity_centroid_mm(*moving) - intensity_centroid_mm(*fixed)."
        )
    pairs = _usable_levels(fixed[0].shape, shrink_factors, smoothing_sigmas_vox)
    shrink_factors, smoothing_sigmas_vox = (tuple(v) for v in zip(*pairs, strict=True))

    tx = s.Euler3DTransform() if with_rotation else s.TranslationTransform(3)
    if with_rotation:
        tx.SetComputeZYX(True)  # R = Rz Ry Rx, matching RigidMotion
        tx.SetCenter(tuple(float(c) for c in center_mm))
    if initial_translation_mm is not None:
        t0 = tuple(float(v) for v in initial_translation_mm)
        tx.SetTranslation(t0) if with_rotation else tx.SetOffset(t0)

    reg = s.ImageRegistrationMethod()
    if metric == "mse":
        reg.SetMetricAsMeanSquares()
        reg.SetOptimizerAsRegularStepGradientDescent(
            learningRate=1.0, minStep=1e-4, numberOfIterations=300, relaxationFactor=0.6
        )
    else:
        reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
        reg.SetMetricSamplingStrategy(reg.RANDOM)
        reg.SetMetricSamplingPercentage(sampling_percentage, int(seed))
        reg.SetOptimizerAsRegularStepGradientDescent(
            learningRate=2.0, minStep=1e-3, numberOfIterations=400, relaxationFactor=0.6
        )
    reg.SetInterpolator(s.sitkLinear)
    reg.SetOptimizerScalesFromPhysicalShift()
    reg.SetShrinkFactorsPerLevel(list(shrink_factors))
    reg.SetSmoothingSigmasPerLevel(list(smoothing_sigmas_vox))
    reg.SetSmoothingSigmasAreSpecifiedInPhysicalUnits(False)
    reg.SetInitialTransform(tx, inPlace=True)
    levels: list[tuple[int, str]] = []
    reg.AddCommand(
        s.sitkMultiResolutionIterationEvent,
        lambda: levels.append(
            (int(reg.GetOptimizerIteration()), reg.GetOptimizerStopConditionDescription())
        ),
    )
    reg.Execute(f_img, m_img)
    final_metric = float(reg.GetMetricValue())
    if not np.isfinite(final_metric) or abs(final_metric) > 1e300:
        raise RuntimeError(
            "registration failed: the metric was never evaluated on any valid sample point "
            "(no overlap between the images at the optimizer's transform)"
        )
    final = (int(reg.GetOptimizerIteration()), reg.GetOptimizerStopConditionDescription())
    all_levels = (*levels[1:], final) if len(levels) > 1 else (final,)
    converged = not any("Maximum number of iterations" in text for _, text in all_levels)
    if not converged:
        warnings.warn(
            "registration hit the iteration limit on at least one resolution level; the result "
            "may be poorly converged (see RegistrationResult.levels)",
            stacklevel=2,
        )

    if with_rotation:
        R = np.asarray(tx.GetMatrix(), float).reshape(3, 3)
        motion = RigidMotion.from_matrix(R, tx.GetTranslation(), center_mm)
    else:
        motion = RigidMotion(tuple(tx.GetOffset()), (0.0, 0.0, 0.0), tuple(center_mm))
    return RegistrationResult(
        motion, final_metric, final[0], final[1], converged, overlap, tuple(all_levels)
    )
