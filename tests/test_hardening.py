"""Malformed or suspicious input must fail with a clear message (or warn), never return a
plausible-looking wrong answer. Each case below was found by an independent review."""

import warnings

import numpy as np
import pydicom
import pytest
import yaml

from doseshift.cli import main
from doseshift.dicomrt import (
    discover_cases,
    load_image_series,
    load_rtdose,
    load_rtstruct_masks,
    write_case,
)
from doseshift.dicomrt.reader import rasterize_contours
from doseshift.geometry import Grid
from doseshift.state import DoseGrid, MaskStructure, PlanState

ROLES = {"PTV": "ptv", "CTV": "ctv", "Rectum": "oar", "Bladder": "oar"}


@pytest.fixture
def paths(tmp_path, case):
    return write_case(case, tmp_path / "case")


def _edit(path, fn):
    ds = pydicom.dcmread(path)
    fn(ds)
    ds.save_as(path, enforce_file_format=True)


# ---- reader: silent wrong answers -------------------------------------------------------------
@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
def test_bad_dose_grid_scaling_is_rejected(paths, bad):
    _edit(paths["dose"], lambda ds: setattr(ds, "DoseGridScaling", bad))
    with pytest.raises(ValueError, match="DoseGridScaling"):
        load_rtdose(paths["dose"])


@pytest.mark.parametrize("bad", [0.0, -5.0, float("nan"), float("inf")])
def test_bad_prescription_is_rejected(case, bad):
    with pytest.raises(ValueError, match="prescription_gy"):
        PlanState.from_masks(case.dose, list(case.masks.values()), bad)


def test_duplicate_roi_names_are_rejected(paths, case):
    def dup(ds):
        ds.StructureSetROISequence[1].ROIName = ds.StructureSetROISequence[0].ROIName

    _edit(paths["rtstruct"], dup)
    with pytest.raises(ValueError, match="more than one ROI is named"):
        load_rtstruct_masks(paths["rtstruct"], case.ct_grid, ["PTV", "CTV"])


def test_nonplanar_contours_warn_instead_of_vanishing(paths, case):
    def point(ds):
        ds.ROIContourSequence[0].ContourSequence[0].ContourGeometricType = "POINT"

    _edit(paths["rtstruct"], point)
    with pytest.warns(UserWarning, match="non-planar"):
        load_rtstruct_masks(paths["rtstruct"], case.ct_grid, ["PTV"])


def test_empty_structure_is_rejected_by_name(case):
    empty = MaskStructure("Ghost", "oar", np.zeros_like(case.masks["PTV"].mask), case.ct_grid)
    with pytest.raises(ValueError, match="'Ghost' is empty"):
        PlanState.from_masks(case.dose, [case.masks["PTV"], empty], 60.0)


def test_structure_outside_the_dose_grid_warns(case):
    small = DoseGrid(Grid.centered((10, 10, 10), (3.0, 3.0, 3.0)), np.ones((10, 10, 10)))
    with pytest.warns(UserWarning, match="outside the dose grid"):
        PlanState.from_masks(small, [case.masks["PTV"]], 60.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # the matching phantom raises no warning
        PlanState.from_masks(case.dose, list(case.masks.values()), 60.0)


def test_negative_or_nonfinite_dose_is_rejected(paths):
    # DoseGridScaling > 0 and unsigned pixels cannot be negative, so corrupt the scaling's effect
    _edit(paths["dose"], lambda ds: setattr(ds, "PixelRepresentation", 1))
    # signed reinterpretation of large unsigned values yields negative doses
    with pytest.raises(ValueError, match="negative|non-finite"):
        load_rtdose(paths["dose"])


# ---- reader: clear errors instead of raw tracebacks ---------------------------------------------
def test_nan_contour_names_the_problem(case):
    poly = np.array([[0.0, 0.0, 0.0], [1.0, np.nan, 0.0], [1.0, 1.0, 0.0]])
    with pytest.raises(ValueError, match="non-finite"):
        rasterize_contours([poly], case.ct_grid, "PTV")


def test_malformed_contour_data_length(paths, case):
    _edit(
        paths["rtstruct"],
        lambda ds: setattr(
            ds.ROIContourSequence[0].ContourSequence[0], "ContourData", [0, 0, 0, 1, 1]
        ),
    )
    with pytest.raises(ValueError, match="malformed ContourData"):
        load_rtstruct_masks(paths["rtstruct"], case.ct_grid, ["PTV"])


def test_unknown_roi_reference(paths, case):
    _edit(
        paths["rtstruct"], lambda ds: setattr(ds.ROIContourSequence[0], "ReferencedROINumber", 9999)
    )
    with pytest.raises(ValueError, match="unknown ROI number 9999"):
        load_rtstruct_masks(paths["rtstruct"], case.ct_grid)


@pytest.mark.parametrize("tag", ["ImageOrientationPatient", "PixelSpacing", "DoseGridScaling"])
def test_missing_required_tags_name_the_file_and_tag(paths, tag):
    _edit(paths["dose"], lambda ds: delattr(ds, tag))
    with pytest.raises(ValueError, match=tag):
        load_rtdose(paths["dose"])


def test_slices_with_different_size_are_rejected(paths):
    f = sorted((paths["ct_dir"]).glob("*.dcm"))[3]
    ds = pydicom.dcmread(f)
    ds.Rows = int(ds.Rows) - 1
    ds.save_as(f, enforce_file_format=True)
    with pytest.raises(ValueError, match="different Rows/Columns"):
        load_image_series(paths["ct_dir"])


def test_missing_directory_and_file_messages(tmp_path):
    with pytest.raises(FileNotFoundError, match="image directory not found"):
        load_image_series(tmp_path / "nope")
    with pytest.raises(FileNotFoundError, match="RTDOSE file not found"):
        load_rtdose(tmp_path / "nope.dcm")


def test_dose_offset_vector_origin_follows_the_standard(paths, case):
    """Frame positions are IPP + normal * GridFrameOffsetVector[k], including the first entry."""

    def shift(ds):
        ds.GridFrameOffsetVector = [float(o) + 6.0 for o in ds.GridFrameOffsetVector]

    ref, _ = load_rtdose(paths["dose"])
    _edit(paths["dose"], shift)
    moved, _ = load_rtdose(paths["dose"])
    assert moved.grid.origin_xyz[2] - ref.grid.origin_xyz[2] == pytest.approx(6.0)


# ---- discovery ----
def test_discover_skips_non_dicom_and_finds_extensionless_files(paths, tmp_path):
    root = tmp_path / "case"
    (root / "notes.dcm").write_text("not dicom")  # first *.dcm in sort order is garbage
    for f in root.rglob("*.dcm"):
        if f.name != "notes.dcm":
            f.rename(f.with_suffix(""))  # TPS exports often have no extension
    (c,) = discover_cases(root)
    assert c.complete and c.prescription_gy == 60.0


def test_discover_missing_directory(tmp_path):
    with pytest.raises(FileNotFoundError, match="directory not found"):
        discover_cases(tmp_path / "nowhere")


# ---- evaluation / experiments ----
def test_empty_cloud_error_names_the_structure(state):
    from doseshift.evaluate import evaluate_doses
    from doseshift.state import PointCloud

    empty = {n: PointCloud(np.empty((0, 3)), np.empty(0)) for n in state.names}
    broken = state.evolve(clouds=empty)
    with pytest.raises(ValueError, match="has no sample points"):
        evaluate_doses(broken, {n: np.empty(0) for n in state.names})


def test_contour_variation_equivalent_is_nan_when_volume_changed(case, state):
    from scipy import ndimage

    from doseshift.experiments import contour_variation

    s = case.masks["Rectum"]
    grown = MaskStructure("Rectum", s.role, ndimage.binary_dilation(s.mask, iterations=3), s.grid)
    df = contour_variation(state, {"grown": {"Rectum": grown}}, metrics={"Rectum": ["Dmean"]})
    assert np.isnan(df["equivalent_shift_mm"]).all() and not df["rigid_comparable"].any()


def test_dimensionless_metrics_use_a_smaller_noise_floor():
    from doseshift.experiments import noise_floor

    assert noise_floor("D95%", 0.1) == 0.1 and noise_floor("V20Gy", 0.1) == 0.1
    assert noise_floor("HI", 0.1) == pytest.approx(0.005)


# ---- CLI ----
def test_cli_reports_errors_without_tracebacks(tmp_path, capsys, paths):
    cases = [
        (["run", str(tmp_path / "nonexistent.yaml")], "config not found"),
        (["demo", "--out", str(tmp_path / "d"), "--courses", "0"], None),  # argparse error
        (["discover", str(tmp_path / "nowhere")], "directory not found"),
        (["discover", str(paths["ct_dir"].parent), "--roles", "PTV"], "--write-config"),
        (
            ["discover", str(paths["ct_dir"].parent), "--roles", "PTV", "--write-config", "x.yaml"],
            "bad --roles",
        ),
        (
            [
                "discover",
                str(paths["ct_dir"].parent),
                "--roles",
                "PTV=ptv",
                "--write-config",
                str(tmp_path / "x.yaml"),
                "--case",
                "5",
            ],
            "out of range",
        ),
    ]
    for argv, needle in cases:
        if needle is None:
            with pytest.raises(SystemExit):
                main(argv)
            continue
        assert main(argv) == 2
        assert needle in capsys.readouterr().err


def test_cli_config_with_null_case_and_bad_n_courses(tmp_path, capsys):
    for body, needle in (
        ({"experiment": "sweep", "case": None}, "case.type"),
        (
            {
                "experiment": "monte_carlo",
                "case": {"type": "synthetic"},
                "monte_carlo": {"n_courses": 0},
            },
            "n_courses",
        ),
    ):
        f = tmp_path / "c.yaml"
        f.write_text(yaml.safe_dump(body))
        assert main(["run", str(f), "--out", str(tmp_path / "o")]) == 2
        assert needle in capsys.readouterr().err


def test_cli_relative_data_paths_give_a_hint(tmp_path, capsys):
    cfg = {
        "experiment": "sweep",
        "case": {
            "type": "dicom",
            "image_dir": "data/raw/missing",
            "dose": "d.dcm",
            "rtstruct": "r.dcm",
            "prescription_gy": 13.0,
            "roles": {"TV": "ptv"},
        },
    }
    f = tmp_path / "c.yaml"
    f.write_text(yaml.safe_dump(cfg))
    assert main(["run", str(f), "--out", str(tmp_path / "o")]) == 2
    assert "run from the repository root" in capsys.readouterr().err


def test_cli_version_and_disclaimer(capsys):
    with pytest.raises(SystemExit):
        main(["--version"])
    assert "doseshift" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        main(["--help"])
    assert "Not validated for clinical use" in capsys.readouterr().out


def test_provenance_rerun_and_dicom_inputs_and_no_anatomy_by_default(tmp_path, paths):
    from doseshift.runner import run_experiment

    cfg = {
        "experiment": "sweep",
        "case": {
            "type": "dicom",
            "image_dir": str(paths["ct_dir"]),
            "dose": str(paths["dose"]),
            "rtstruct": str(paths["rtstruct"]),
            "prescription_gy": 60.0,
            "roles": ROLES,
        },
        "sweep": {"max_mm": 1, "step_mm": 1},
    }
    files = run_experiment(cfg, tmp_path / "run1")
    import json

    prov = json.loads(files["provenance.json"].read_text())
    assert prov["resolved"]["show_image"] is False  # real data: no anatomy unless requested
    assert {"matplotlib", "pydicom", "numpy"} <= set(prov["packages"]) and "platform" in prov
    assert len(prov["inputs"]["dose"]["sha256"]) == 64 and prov["inputs"]["image_slices"] == 48
    # a provenance file is itself a valid config: it re-runs the same experiment
    assert main(["run", str(files["provenance.json"]), "--out", str(tmp_path / "run2")]) == 0
    a = (tmp_path / "run1" / "sweep.csv").read_text()
    assert a == (tmp_path / "run2" / "sweep.csv").read_text()


# ---- findings of the scientific review ----
def test_structure_sampling_conserves_volume_when_part_of_the_structure_is_outside_the_dose_grid(
    case,
):
    """Sampling a structure on the dose grid must not redefine it: the part outside the grid stays
    (reading 0 Gy, like native evaluation) instead of being silently dropped."""
    from doseshift.perturb import StructureSampling

    small = DoseGrid(Grid.centered((12, 12, 12), (3.0, 3.0, 3.0)), np.full((12, 12, 12), 60.0))
    with pytest.warns(UserWarning, match="outside the dose grid"):
        st = PlanState.from_masks(small, [case.masks["PTV"]], 60.0)
    v0 = st.clouds["PTV"].total_volume
    for mode in ("dose_center", "dose_fractional"):
        v = StructureSampling(mode)(st).clouds["PTV"].total_volume
        assert abs(v - v0) / v0 < 0.06, mode


def test_structure_sampling_margin_covers_structures_with_coarse_voxels_on_a_fine_dose_grid(case):
    from doseshift.perturb import StructureSampling

    fine = DoseGrid(Grid.centered((160, 160, 160), (0.5, 0.5, 0.5)), np.full((160,) * 3, 60.0))
    st = PlanState.from_masks(fine, [case.masks["CTV"]], 60.0)
    v0 = st.clouds["CTV"].total_volume
    v = StructureSampling("dose_center")(st).clouds["CTV"].total_volume
    assert abs(v - v0) / v0 < 0.01  # before the fix a too-small candidate box lost ~1 %


def test_coarsen_warns_when_it_drops_trailing_voxels(state):
    from doseshift.perturb import DoseCoarsen

    odd = DoseGrid(Grid.centered((47, 61, 61), (3.0, 3.0, 3.0)), np.ones((47, 61, 61)))
    st = PlanState.from_masks(odd, [state.masks["PTV"]], 60.0)
    with pytest.warns(UserWarning, match="not a multiple"):
        DoseCoarsen(2)(st)


def test_dose_offset_equals_patient_shift_minus_A_v_on_oblique_negative_k_grids(case):
    """The documented relation is t = -A v with A the affine's 3x3 part, not -v * spacing."""
    from scipy.spatial.transform import Rotation

    from doseshift.evaluate import evaluate
    from doseshift.perturb import DoseGridOffset, RigidShift

    R = Rotation.from_euler("z", 20, degrees=True).as_matrix()
    direction = R @ np.diag([1.0, 1.0, -1.0])  # oblique in-plane, slice direction reversed
    shape, spacing = (40, 50, 50), (2.0, 2.5, 3.0)  # (nz, ny, nx), (sx, sy, sz)
    A = direction @ np.diag(spacing)
    centre_idx = np.array([(shape[2] - 1) / 2, (shape[1] - 1) / 2, (shape[0] - 1) / 2])
    origin = np.array([0.0, 40.0, 0.0]) - A @ centre_idx  # grid centred on the rectum
    grid = Grid.from_spacing(shape, spacing, tuple(origin), direction)
    kk, jj, ii = np.meshgrid(*(np.arange(n) for n in shape), indexing="ij")
    xyz = grid.index_to_world(np.stack([ii.ravel(), jj.ravel(), kk.ravel()], axis=1))
    dose = 60.0 * np.exp(-(((xyz - [0.0, 20.0, 0.0]) / 30.0) ** 2).sum(axis=1) / 2)
    st = PlanState.from_masks(DoseGrid(grid, dose.reshape(shape)), [case.masks["Rectum"]], 60.0)
    assert st.outside_dose_fraction("Rectum") < 0.01  # the comparison below must not be vacuous

    v = np.array([0.5, -0.25, 1.0])
    t = -(A @ v)
    keys = {"Rectum": ["Dmean", "D2%"]}
    a = evaluate(RigidShift(tuple(t))(st), keys).flat()
    b = evaluate(DoseGridOffset(tuple(v))(st), keys).flat()
    assert b["Rectum:Dmean"] > 1.0  # real dose there
    for k in a:
        assert a[k] == pytest.approx(b[k], rel=1e-9, abs=1e-9)
    naive = evaluate(RigidShift(tuple(-v * grid.spacing_xyz))(st), {"Rectum": ["Dmean"]}).flat()
    assert abs(naive["Rectum:Dmean"] - b["Rectum:Dmean"]) > 1e-3  # the old relation is wrong here
