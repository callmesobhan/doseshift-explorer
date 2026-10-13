import numpy as np

from doseshift.geometry import Grid, RigidMotion


def test_index_world_roundtrip_oblique():
    d = np.array([[0.98, -0.199, 0.0], [0.199, 0.98, 0.0], [0.0, 0.0, 1.0]])
    g = Grid.from_spacing((5, 6, 7), (1.5, 2.0, 3.0), (10, -20, 5), direction=d)
    ijk = np.array([[0, 0, 0], [3.5, 2.0, 1.25], [6, 5, 4]], float)
    assert np.allclose(g.world_to_index(g.index_to_world(ijk)), ijk)
    assert np.allclose(g.spacing_xyz, (1.5, 2.0, 3.0))
    assert np.isclose(g.voxel_volume, 1.5 * 2.0 * 3.0)


def test_centered_grid_is_symmetric():
    g = Grid.centered((4, 6, 8), (2.0, 3.0, 5.0))
    lo, hi = g.extent_xyz()
    assert np.allclose(lo, -hi)


def test_mask_points_match_voxel_centres():
    g = Grid.from_spacing((3, 3, 3), (2.0, 2.0, 2.0), (100, 200, 300))
    m = np.zeros(g.shape, bool)
    m[1, 2, 0] = True  # k=1, j=2, i=0
    assert np.allclose(g.mask_points(m), [[100.0, 204.0, 302.0]])


def test_rotation_fixes_center_and_is_right_handed():
    m = RigidMotion(rotation_deg=(0, 0, 90), center_mm=(10, 0, 0))
    out = m.apply(np.array([[10.0, 0.0, 0.0], [11.0, 0.0, 0.0]]))
    assert np.allclose(out[0], [10, 0, 0])
    assert np.allclose(out[1], [10, 1, 0])  # +x rotates onto +y for +90 deg about z


def test_inverse_apply_undoes_apply():
    m = RigidMotion((1, -2, 3), (4, 5, -6), (7, 8, 9))
    p = np.random.default_rng(0).normal(size=(20, 3)) * 30
    assert np.allclose(m.inverse_apply(m.apply(p)), p)


def test_warp_volume_moves_features_by_translation():
    g = Grid.from_spacing((20, 20, 20), (1.0, 1.0, 1.0))
    vol = np.zeros(g.shape)
    vol[10, 10, 8] = 1.0
    out = RigidMotion((3.0, 0.0, 0.0)).warp_volume(g, vol, order=1)
    assert np.unravel_index(out.argmax(), out.shape) == (10, 10, 11)
