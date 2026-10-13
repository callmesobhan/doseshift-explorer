"""Validate the motion convention against an independent implementation (SimpleITK)."""

import numpy as np
import pytest

pytest.importorskip("SimpleITK")

from doseshift.geometry import RigidMotion  # noqa: E402
from doseshift.registration import register_rigid  # noqa: E402
from helpers import cbct_like  # noqa: E402


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


# ---- cross-scanner: different grid, frame origin and intensity mapping ----
@pytest.mark.parametrize(
    "t,r,start",
    [
        ((3.0, -2.0, 4.0), (0, 0, 0), (0.0, 0.0, 0.0)),
        ((-4.0, 3.0, 2.0), (1.0, -1.5, 2.0), (12.0, 0.0, 0.0)),
        ((5.0, 5.0, -3.0), (0.0, 2.0, -2.5), (0.0, -12.0, 0.0)),
    ],
)
def test_mutual_information_recovers_motion_across_grid_frame_and_intensity(case, t, r, start):
    from doseshift.registration import intensity_centroid_mm

    fg, fv = case.ct_grid, case.ct_hu
    lo, hi = fg.extent_xyz()
    centre = tuple((lo + hi) / 2)
    offset = np.array([83.0, -41.0, 117.0])
    mg, mv = cbct_like(fg, fv, RigidMotion(t, r, centre), offset, seed=1)
    t0 = intensity_centroid_mm(mg, mv) - intensity_centroid_mm(fg, fv) + np.array(start)
    res = register_rigid(
        (fg, fv), (mg, mv), center_mm=centre, metric="mi", initial_translation_mm=tuple(t0)
    )
    # the frame offset is part of the translation: anatomy at p sits at R(p-c)+c+t+offset
    assert np.allclose(res.motion.translation_mm, np.array(t) + offset, atol=0.5)
    assert np.allclose(res.motion.rotation_deg, r, atol=1.0)


def test_centroid_start_is_independent_of_the_intensity_scale(case):
    from doseshift.registration import intensity_centroid_mm, otsu_threshold

    fg, fv = case.ct_grid, case.ct_hu
    ref = intensity_centroid_mm(fg, fv)
    assert np.allclose(intensity_centroid_mm(fg, 0.37 * fv + 220.0), ref, atol=1e-6)
    assert np.allclose(ref, [0.0, 0.0, 0.0], atol=2.0)  # phantom body is centred on the origin
    assert -1000 < otsu_threshold(fv) < 200


def test_register_rigid_rejects_unknown_metric(case):
    with pytest.raises(ValueError, match="metric"):
        register_rigid((case.ct_grid, case.ct_hu), (case.ct_grid, case.ct_hu), metric="ssd")


# ---- failures must be loud, never an identity "perfect setup" (independent review) ----
def test_unrelated_frames_raise_instead_of_returning_identity(case):
    from doseshift.geometry import Grid

    fg, fv = case.ct_grid, case.ct_hu
    far = Grid(fg.shape, fg.affine.copy())
    far.affine[:3, 3] += [1000.0, 0.0, 0.0]  # same image, coordinates 1 m away
    with pytest.raises(ValueError, match="barely overlap"):
        register_rigid((fg, fv), (far, fv))


def test_non_finite_images_are_rejected(case):
    bad = case.ct_hu.copy()
    bad[3, 3, 3] = np.nan
    with pytest.raises(ValueError, match="NaN or infinite"):
        register_rigid((case.ct_grid, case.ct_hu), (case.ct_grid, bad))


def test_result_reports_convergence_and_every_level(case):
    truth = RigidMotion((3.0, 0.0, 0.0), (0, 0, 0))
    moving = truth.warp_volume(case.ct_grid, case.ct_hu, cval=-1000.0)
    res = register_rigid((case.ct_grid, case.ct_hu), (case.ct_grid, moving))
    assert res.converged and res.overlap_fraction > 0.9
    assert len(res.levels) == 3 and res.iterations == res.levels[-1][0]


def test_very_thin_volumes_get_a_clear_error_not_a_raw_itk_exception(case):
    thin_grid = type(case.ct_grid)((3, *case.ct_grid.shape[1:]), case.ct_grid.affine)
    thin = np.ascontiguousarray(case.ct_hu[20:23])
    with pytest.raises(ValueError, match="at least 4"):
        register_rigid((thin_grid, thin), (thin_grid, thin))


# ---- independent check: the moving image comes from SimpleITK's own resampler -----------------
def _sitk_moved(grid, values, motion, default=-1000.0):
    """The patient after ``motion``, made by SimpleITK (not by RigidMotion.warp_volume / scipy)."""
    import SimpleITK as sitk

    from doseshift.registration import to_sitk

    img = to_sitk(grid, values)
    R = motion.matrix
    inv = sitk.AffineTransform(3)  # x -> Rinv (x - c) + c - Rinv t  == the inverse of the motion
    inv.SetMatrix(tuple(R.T.ravel()))
    inv.SetCenter(tuple(float(c) for c in motion.center_mm))
    inv.SetTranslation(tuple(float(v) for v in -(R.T @ np.asarray(motion.translation_mm))))
    out = sitk.Resample(img, img, inv, sitk.sitkLinear, default)
    return sitk.GetArrayFromImage(out)


def _unusual_grid(kind):
    from scipy.spatial.transform import Rotation

    from doseshift.geometry import Grid

    shape, sp = (54, 62, 62), (3.0, 3.0, 3.0)  # +-80..92 mm: the whole phantom body fits
    rot = Rotation.from_euler("zx", [20, 12], degrees=True).as_matrix()
    flip_k, flip_i = np.diag([1.0, 1.0, -1.0]), np.diag([-1.0, 1.0, 1.0])
    d = {
        "oblique": rot,
        "negative-k": flip_k,
        "left-handed": flip_i,
        "oblique+negative-k": rot @ flip_k,
    }[kind]
    A = d @ np.diag(sp)
    centre_idx = np.array([(shape[2] - 1) / 2, (shape[1] - 1) / 2, (shape[0] - 1) / 2])
    return Grid.from_spacing(shape, sp, tuple(-A @ centre_idx), d)


@pytest.mark.parametrize("kind", ["oblique", "negative-k", "left-handed", "oblique+negative-k"])
def test_registration_agrees_with_an_independent_resampler_on_unusual_grids(kind):
    from doseshift.synthetic import make_ct

    grid = _unusual_grid(kind)
    values = make_ct(grid)
    centre = (5.0, -4.0, 3.0)  # off-centre rotation point
    truth = RigidMotion((3.0, -2.0, 2.5), (1.5, -1.0, 2.0), centre)
    moving = _sitk_moved(grid, values, truth)
    res = register_rigid((grid, values), (grid, moving), center_mm=centre)
    pts = np.random.default_rng(0).uniform(-35, 35, (200, 3))
    assert np.abs(res.motion.apply(pts) - truth.apply(pts)).max() < 0.3


def test_registration_degrades_gracefully_when_the_field_of_view_truncates_the_body(case):
    """The accuracy quoted for full-field images is a best case; real CT/CBCT are truncated."""
    from doseshift.geometry import Grid

    g, v = case.ct_grid, case.ct_hu
    crop = (slice(None), slice(21, 69), slice(21, 69))  # +-47 mm: cuts through the body
    cg = Grid(v[crop].shape, g.translated_voxels((21, 21, 0)).affine)
    cv = np.ascontiguousarray(v[crop])
    truth = RigidMotion((3.0, -2.0, 2.0), (1.0, 0.0, 1.5), (0.0, 0.0, 0.0))
    res = register_rigid((cg, cv), (cg, _sitk_moved(cg, cv, truth)), center_mm=(0.0, 0.0, 0.0))
    pts = np.random.default_rng(1).uniform(-30, 30, (200, 3))
    assert np.abs(res.motion.apply(pts) - truth.apply(pts)).max() < 1.5  # worse than full field
