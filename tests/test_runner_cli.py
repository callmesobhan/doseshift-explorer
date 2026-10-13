from pathlib import Path

import pandas as pd
import pytest

from doseshift.cli import main
from doseshift.config import EXPERIMENTS, load_config, validate_config
from doseshift.runner import demo_configs, run_experiment

EXPERIMENT_DIR = Path(__file__).parent.parent / "experiments"


@pytest.mark.parametrize("path", sorted(EXPERIMENT_DIR.glob("*.yaml")), ids=lambda p: p.name)
def test_example_configs_parse(path):
    cfg = load_config(path)
    assert cfg["experiment"] in EXPERIMENTS


def test_validation_errors():
    with pytest.raises(ValueError):
        validate_config({"experiment": "nope", "case": {"type": "synthetic"}})
    with pytest.raises(ValueError):
        validate_config({"experiment": "sweep", "case": {"type": "dicom"}})
    with pytest.raises(ValueError):
        validate_config(
            {
                "experiment": "margin_scan",
                "case": {
                    "type": "dicom",
                    "ct_dir": 1,
                    "dose": 1,
                    "rtstruct": 1,
                    "prescription_gy": 1,
                    "roles": {},
                },
            }
        )


def test_each_experiment_writes_tables_figures_and_provenance(tmp_path):
    cfgs = demo_configs(n_courses=6)
    cfgs["sweep"]["sweep"].update(max_mm=2.0, step_mm=1.0)
    cfgs["monte_carlo"]["monte_carlo"]["error_model"]["n_fractions"] = 2
    cfgs["margin_scan"]["margin_scan"].update(margins_mm=[4, 7], n_courses=20)
    for name, cfg in cfgs.items():
        files = run_experiment(cfg, tmp_path / name)
        assert files["provenance.json"].exists() and files["geometry.png"].stat().st_size > 1000
        for p in files.values():
            assert p.exists() and p.stat().st_size > 0, (name, p)
    assert "equivalent_shift_mm" in pd.read_csv(tmp_path / "grid_audit" / "audit_equivalence.csv")


def test_cli_phantom_inspect_and_dicom_run(tmp_path, capsys):
    dcm = tmp_path / "phantom"
    assert main(["make-phantom", str(dcm)]) == 0
    assert (
        main(
            [
                "inspect",
                "--ct-dir",
                str(dcm / "CT"),
                "--dose",
                str(dcm / "RTDOSE.dcm"),
                "--rtstruct",
                str(dcm / "RTSTRUCT.dcm"),
            ]
        )
        == 0
    )
    assert "CT == DOSE : True" in capsys.readouterr().out
    cfg = {
        "experiment": "sweep",
        "case": {
            "type": "dicom",
            "ct_dir": str(dcm / "CT"),
            "dose": str(dcm / "RTDOSE.dcm"),
            "rtstruct": str(dcm / "RTSTRUCT.dcm"),
            "prescription_gy": 60.0,
            "roles": {"PTV": "ptv", "CTV": "ctv", "Rectum": "oar", "Bladder": "oar"},
        },
        "sweep": {"max_mm": 2, "step_mm": 1},
    }
    files = run_experiment(cfg, tmp_path / "dicom_sweep")
    assert files["sweep.csv"].exists() and files["sweep.png"].exists()
