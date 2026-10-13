"""Rigid image registration (SimpleITK) expressed in the project's motion convention.

Registering a *treatment* image (moving) to the *planning* CT (fixed) yields a transform that
maps planned positions to where the same anatomy sits on the day: exactly ``RigidMotion``. That
lets measured setup errors feed the Monte Carlo error model directly, and lets the perturbation
engine be validated against an independent implementation (tests/test_registration.py).

Requires the ``registration`` extra: ``pip install doseshift[registration]``.
"""

from __future__ import annotations

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
    motion: RigidMotion
    metric_value: float
    iterations: int
    stop_condition: str


def register_rigid(
    fixed: tuple[Grid, np.ndarray],
    moving: tuple[Grid, np.ndarray],
    center_mm: tuple[float, float, float] | None = None,
    with_rotation: bool = True,
    shrink_factors: tuple[int, ...] = (4, 2, 1),
    smoothing_sigmas_vox: tuple[float, ...] = (2.0, 1.0, 0.0),
) -> RegistrationResult:
    """Mean-squares rigid registration (same-modality CT/CBCT-like images).

    ``center_mm`` is the rotation centre (use the PTV centroid so rotation and translation are
    reported about the isocentre); default is the fixed image centre.
    """
    s = _require_sitk()
    f_img, m_img = to_sitk(*fixed), to_sitk(*moving)
    if center_mm is None:
        lo, hi = fixed[0].extent_xyz()
        center_mm = tuple((lo + hi) / 2)

    tx = s.Euler3DTransform() if with_rotation else s.TranslationTransform(3)
    if with_rotation:
        tx.SetComputeZYX(True)  # R = Rz Ry Rx, matching RigidMotion
        tx.SetCenter(tuple(float(c) for c in center_mm))

    reg = s.ImageRegistrationMethod()
    reg.SetMetricAsMeanSquares()
    reg.SetInterpolator(s.sitkLinear)
    reg.SetOptimizerAsRegularStepGradientDescent(
        learningRate=1.0, minStep=1e-4, numberOfIterations=300, relaxationFactor=0.6
    )
    reg.SetOptimizerScalesFromPhysicalShift()
    reg.SetShrinkFactorsPerLevel(list(shrink_factors))
    reg.SetSmoothingSigmasPerLevel(list(smoothing_sigmas_vox))
    reg.SetSmoothingSigmasAreSpecifiedInPhysicalUnits(False)
    reg.SetInitialTransform(tx, inPlace=True)
    reg.Execute(f_img, m_img)

    if with_rotation:
        R = np.asarray(tx.GetMatrix(), float).reshape(3, 3)
        motion = RigidMotion.from_matrix(R, tx.GetTranslation(), center_mm)
    else:
        motion = RigidMotion(tuple(tx.GetOffset()), (0.0, 0.0, 0.0), tuple(center_mm))
    return RegistrationResult(
        motion,
        float(reg.GetMetricValue()),
        int(reg.GetOptimizerIteration()),
        reg.GetOptimizerStopConditionDescription(),
    )
