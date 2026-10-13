import shutil
import warnings

import numpy as np
import pydicom
import pytest
import yaml

from doseshift.cli import main
from doseshift.dicomrt import (
    case_config,
    discover_cases,
    load_image_series,
    load_plan_state,
    write_case,
)
from doseshift.evaluate import evaluate
from doseshift.synthetic import make_case

ROLES = {"PTV": "ptv", "CTV": "ctv", "Rectum": "oar", "Bladder": "oar"}


@pytest.fixture(scope="module")
def flat(tmp_path_factory, case):
    root = tmp_path_factory.mktemp("flat")
    return root, write_case(case, root)


def test_discover_follows_dicom_references_in_a_flat_folder(flat):
    root, paths = flat
    (c,) = discover_cases(root)
    assert c.image_dir == paths["ct_dir"]
    assert c.rtstruct == paths["rtstruct"] and c.rtdose == paths["dose"]
    assert c.rtplan == paths["rtplan"]
    assert c.prescription_gy == 60.0 and c.n_fractions == 20
    assert set(c.roi_names) == set(ROLES) and c.complete and not c.notes


def test_discover_without_rtplan_pairs_by_frame_of_reference_and_has_no_prescription(
    flat, tmp_path
):
    root, _ = flat
    shutil.copytree(root, tmp_path / "noplan")
    (tmp_path / "noplan" / "RTPLAN.dcm").unlink()
    (c,) = discover_cases(tmp_path / "noplan")
    assert c.rtdose is not None and c.rtplan is None and c.prescription_gy is None
    assert any("FrameOfReference" in n for n in c.notes)
    with pytest.raises(ValueError, match="prescription unknown"):
        case_config(c, ROLES)
    cfg = case_config(c, ROLES, prescription_gy=50.0)
    assert cfg["case"]["prescription_is_reference"] is True
    assert "CI" not in cfg["metrics"]["PTV"]  # Rx-dependent metrics dropped


def test_case_config_runs_end_to_end(flat, tmp_path):
    root, _ = flat
    out = tmp_path / "cfg.yaml"
    assert (
        main(["discover", str(root), "--roles", "PTV=ptv,Rectum=oar", "--write-config", str(out)])
        == 0
    )
    cfg = yaml.safe_load(out.read_text())
    assert cfg["case"]["prescription_gy"] == 60.0 and "prescription_is_reference" not in cfg["case"]
    assert main(["run", str(out), "--out", str(tmp_path / "run")]) == 0


def test_unreferenced_rtstruct_triggers_a_warning(case, flat, tmp_path):
    root, paths = flat
    other = write_case(make_case(seed=0), tmp_path / "other")  # different series UIDs
    with pytest.warns(UserWarning, match="references image series"):
        load_plan_state(other["ct_dir"], paths["dose"], paths["rtstruct"], ROLES, 60.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # matching pair: no warnings at all
        load_plan_state(paths["ct_dir"], paths["dose"], paths["rtstruct"], ROLES, 60.0)


def test_mr_series_load_like_ct(flat, tmp_path):
    _, paths = flat
    mr = tmp_path / "MR"
    shutil.copytree(paths["ct_dir"], mr)
    for f in mr.glob("*.dcm"):
        ds = pydicom.dcmread(f)
        ds.Modality = "MR"
        ds.save_as(f, enforce_file_format=True)
    grid, values, _ = load_image_series(mr)
    assert grid.shape == make_case().ct_grid.shape and np.isfinite(values).all()
    with pytest.raises(FileNotFoundError):
        load_image_series(mr, modalities=("CT",))


def test_reference_dose_label_is_carried(flat):
    _, paths = flat
    state = load_plan_state(
        paths["ct_dir"], paths["dose"], paths["rtstruct"], ROLES, 60.0, rx_label="reference dose"
    )
    assert state.rx_label == "reference dose"
    assert evaluate(state).flat()["PTV:D95%"] > 0
