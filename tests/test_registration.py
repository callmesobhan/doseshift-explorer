"""Validate the motion convention against an independent implementation (SimpleITK)."""

import numpy as np
import pytest

pytest.importorskip("SimpleITK")

from doseshift.geometry import RigidMotion  # noqa: E402
from doseshift.registration import register_rigid  # noqa: E402


@pytest.mark.parametrize(
    "translation,rotation",
    [
        ((3.0, 0.0, 0.0), (0, 0, 0)),
        ((2.5, -1.5, 3.0), (1.0, 0.0, 2.0)),
        ((0, 4.0, -2.0), (0, 2.0, 0)),
    ],
)
def test_registration_recovers_injected_patient_motion(case, state, translation, rotation):
    iso = tuple(state.isocenter())
    truth = RigidMotion(translation, rotation, iso)
    moving = truth.warp_volume(case.ct_grid, case.ct_hu, order=1, cval=-1000.0)

    result = register_rigid((case.ct_grid, case.ct_hu), (case.ct_grid, moving), center_mm=iso)

    assert np.allclose(result.motion.translation_mm, truth.translation_mm, atol=0.4)
    assert np.allclose(result.motion.rotation_deg, truth.rotation_deg, atol=0.4)
    # the recovered motion must act on structure points like the injected one
    pts = state.clouds["Rectum"].points[::50]
    assert np.abs(result.motion.apply(pts) - truth.apply(pts)).max() < 0.8


def test_translation_only_mode(case):
    truth = RigidMotion((2.0, 1.0, -3.0))
    moving = truth.warp_volume(case.ct_grid, case.ct_hu, order=1, cval=-1000.0)
    result = register_rigid((case.ct_grid, case.ct_hu), (case.ct_grid, moving), with_rotation=False)
    assert np.allclose(result.motion.translation_mm, truth.translation_mm, atol=0.4)
