"""Run a configured experiment and write CSV tables, figures and provenance to a directory.

``run_experiment`` returns ``{filename: path}`` for everything it wrote."""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__, viz
from .config import (
    build_state,
    constraints_from,
    error_model_from,
    validate_config,
)
from .evaluate import evaluate
from .experiments import (
    audit,
    audit_equivalence,
    contour_variation,
    directional_sweeps,
    margin_scan,
)
from .montecarlo import ErrorModel, simulate
from .perturb import DoseCoarsen, DoseGridOffset, RigidShift

DEFAULT_KEYS = ["CTV:D98%", "PTV:D95%", "Rectum:Dmean"]


def _versions() -> dict[str, str]:
    import importlib.metadata as md

    out = {}
    for pkg in ("numpy", "scipy", "pandas", "matplotlib", "pydicom", "SimpleITK", "PyYAML"):
        try:
            out[pkg] = md.version(pkg)
        except md.PackageNotFoundError:
            pass
    return out


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _dicom_inputs(case: dict) -> dict:
    """Identify the exact input objects of a DICOM run (UIDs and file hashes, no patient data)."""
    import pydicom

    from .config import image_dir

    out: dict = {}
    img = Path(image_dir(case))
    uid = None
    n = 0
    for f in sorted(img.iterdir()):
        if f.is_file():
            try:
                ds = pydicom.dcmread(f, stop_before_pixels=True)
            except pydicom.errors.InvalidDicomError:
                continue
            uid, n = str(ds.SeriesInstanceUID), n + 1
    out["image_series_uid"], out["image_slices"] = uid, n
    for key in ("dose", "rtstruct"):
        p = Path(case[key])
        out[key] = {
            "sop_instance_uid": str(pydicom.dcmread(p, stop_before_pixels=True).SOPInstanceUID),
            "sha256": _sha256(p),
        }
    return out


def _provenance(cfg: dict, outdir: Path, resolved: dict) -> Path:
    path = outdir / "provenance.json"
    info = {
        "doseshift": __version__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": _versions(),
        "resolved": resolved,
        "config": cfg,
    }
    if cfg["case"]["type"] == "dicom":
        info["inputs"] = _dicom_inputs(cfg["case"])
    path.write_text(json.dumps(info, indent=2, default=str))
    return path


def run_experiment(cfg: dict, outdir: str | Path) -> dict[str, Path]:
    validate_config(cfg)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    theme = cfg.get("theme", "light")
    state, ct = build_state(cfg["case"])
    metrics = cfg.get("metrics")
    kind = cfg["experiment"]
    files: dict[str, Path] = {}

    def table(name: str, df: pd.DataFrame, index: bool = False) -> None:
        files[f"{name}.csv"] = out / f"{name}.csv"
        df.to_csv(files[f"{name}.csv"], index=index)

    def figure(name: str, fig) -> None:
        files[f"{name}.png"] = out / f"{name}.png"
        viz.save(fig, files[f"{name}.png"])

    synthetic = cfg["case"]["type"] == "synthetic"
    title = cfg.get("name") or ("Phantom plan geometry" if synthetic else "Plan geometry")
    # Real image slices are not rendered unless asked for (show_image: true): outputs of public
    # repositories and CI must not embed patient anatomy by accident.
    show_image = bool(cfg.get("show_image", synthetic))
    figure(
        "geometry",
        viz.plot_slices(state, ct=ct if show_image else None, theme=theme, title=title),
    )
    table("nominal", evaluate(state, metrics).to_frame())

    if kind == "sweep":
        sc = cfg.get("sweep", {})
        keys = sc.get("keys", DEFAULT_KEYS)
        sweeps = directional_sweeps(state, sc.get("max_mm", 10.0), sc.get("step_mm", 0.5), metrics)
        table("sweep", pd.concat(sweeps.values(), ignore_index=True))
        figure("sweep", viz.plot_sweeps(sweeps, keys, sc.get("thresholds"), theme))

    elif kind == "grid_audit":
        ac = cfg.get("audit", {})
        keys = ac.get("keys", ["PTV:D95%", "Rectum:D2%"])
        df = audit(state, metrics=metrics)
        eq = audit_equivalence(
            state, df, keys, metrics, ac.get("max_mm", 12.0), ac.get("min_delta", 0.1)
        )
        table("audit", df)
        table("audit_equivalence", eq)
        figure("audit_equivalence", viz.plot_equivalence(eq, keys, ac.get("max_mm", 12.0), theme))
        scenarios = {
            "planned": evaluate(state, metrics),
            "dose offset +1/2 voxel (ijk)": evaluate(DoseGridOffset((0.5,) * 3)(state), metrics),
            "dose grid coarsened x2": evaluate(DoseCoarsen(2)(state), metrics),
            "patient shifted 3 mm (y)": evaluate(RigidShift((0, 3.0, 0))(state), metrics),
        }
        figure(
            "dvh_scenarios",
            viz.plot_dvh_overlay(
                scenarios,
                state.prescription_gy,
                theme=theme,
                rx_label=state.rx_label,
                title="DVH: grid assumptions vs a physical shift",
            ),
        )

    elif kind == "monte_carlo":
        mc = cfg.get("monte_carlo", {})
        model = error_model_from(mc.get("error_model"))
        res = simulate(
            state,
            model,
            mc.get("n_courses", 300),
            cfg.get("seed", 0),
            mc.get("method", "sampled"),
            metrics,
        )
        table("mc_metrics", res.metrics)
        table("mc_summary", res.summary().reset_index())
        table("mc_draws", res.draws)
        figure("mc_dvh_band", viz.plot_mc_bands(res, theme=theme, rx_label=state.rx_label))
        cons = constraints_from(mc.get("constraints"))
        if cons:
            table(
                "mc_probabilities",
                pd.DataFrame(
                    [
                        {
                            "constraint": str(c),
                            "probability": res.probability(c),
                            "se": res.probability_se(c),
                            "n_courses": res.n_courses,
                            "planned_passes": bool(
                                c.passes(
                                    np.array([res.nominal.flat()[c.key]]), state.prescription_gy
                                )[0]
                            ),
                        }
                        for c in cons
                    ]
                ),
            )
            figure("mc_hist", viz.plot_metric_distribution(res, cons[0].key, cons[0], theme))

    elif kind == "margin_scan":
        mc = cfg.get("margin_scan", {})
        model = error_model_from(mc.get("error_model"))
        con = constraints_from([mc["constraint"]])[0]
        from .synthetic import make_state

        seed = int(cfg["case"].get("seed", 0))
        df = margin_scan(
            lambda m: make_state(seed, m),
            mc.get("margins_mm", [2, 3, 4, 5, 6, 7, 8, 9]),
            model,
            con,
            mc.get("n_courses", 500),
            cfg.get("seed", 0),
            mc.get("method", "blur"),
        )
        table("margin_scan", df)
        figure(
            "margin_scan",
            viz.plot_margin_scan(
                df, float(model.van_herk_margin_mm().max()), label=f"P({con})", theme=theme
            ),
        )

    elif kind == "contour_variation":
        _contour_variation(cfg, state, metrics, theme, table, figure)

    resolved = {"seed": cfg.get("seed", 0), "theme": theme, "show_image": show_image}
    files["provenance.json"] = _provenance(cfg, out, resolved)
    return files


def _contour_variation(cfg, state, metrics, theme, table, figure) -> None:
    from .config import image_dir
    from .dicomrt import load_image_series, load_rtstruct_masks
    from .state import MaskStructure, PlanState

    cv = cfg["contour_variation"]
    alternatives = {}
    for alt in cv["alternatives"]:
        grid, _, _ = load_image_series(image_dir(alt))
        masks, _ = load_rtstruct_masks(alt["rtstruct"], grid, list(alt["names"].values()))
        alternatives[alt["label"]] = {
            plan: MaskStructure(plan, state.role(plan), masks[other], grid)
            for plan, other in alt["names"].items()
        }
    names = list(next(iter(alternatives.values())))
    df = contour_variation(
        state,
        alternatives,
        metrics,
        cv.get("max_mm", 12.0),
        cv.get("min_delta", 0.1),
        cv.get("max_volume_change", 0.25),
    )
    table("contour_variation", df)

    rx = state.prescription_gy
    planned = PlanState.from_masks(state.dose, [state.masks[n] for n in names], rx)
    evals = {"planning contours": evaluate(planned, metrics)}
    for label, structs in alternatives.items():
        alt_state = PlanState.from_masks(state.dose, [structs[n] for n in names], rx)
        evals[label] = evaluate(alt_state, metrics)
    figure(
        "dvh_contours",
        viz.plot_dvh_overlay(
            evals,
            rx,
            title="Same planned dose, different contours",
            subtitle="planning contours vs contours drawn on later images",
            theme=theme,
            rx_label=state.rx_label,
        ),
    )

    keys = cv.get("keys") or [f"{n}:D2%" for n in names]
    eq = df.assign(scenario=df["alternative"], metric=df["structure"] + ":" + df["metric"])
    eq = eq[eq["metric"].isin(keys)].copy()
    eq.loc[~eq["rigid_comparable"], "equivalent_shift_mm"] = float("nan")
    figure(
        "contour_equivalence",
        viz.plot_equivalence(
            eq,
            keys,
            cv.get("max_mm", 12.0),
            theme,
            title="How much rigid shift is the measured contour change worth?",
            subtitle="Smallest axis-aligned shift with the same effect (lower bound)",
            na_label="n/a: no rigid shift reproduces this (or organ volume changed)",
        ),
    )


def demo_configs(n_courses: int = 300) -> dict[str, dict]:
    case = {"type": "synthetic", "seed": 0, "ctv_margin_mm": 5.0}
    em = {"sys_sigma_mm": [2, 2, 2], "rand_sigma_mm": [3, 3, 3], "n_fractions": 20}
    con = [{"structure": "CTV", "metric": "D98%", "op": ">=", "value": 95.0, "pct_of_rx": True}]
    return {
        "sweep": {
            "experiment": "sweep",
            "case": case,
            "sweep": {"max_mm": 10, "step_mm": 0.5, "thresholds": {"CTV:D98%": 57.0}},
        },
        "grid_audit": {
            "experiment": "grid_audit",
            "case": case,
            "audit": {"keys": ["PTV:D95%", "Rectum:D2%"]},
        },
        "monte_carlo": {
            "experiment": "monte_carlo",
            "case": case,
            "seed": 1,
            "monte_carlo": {
                "n_courses": n_courses,
                "method": "sampled",
                "error_model": em,
                "constraints": con,
            },
        },
        "margin_scan": {
            "experiment": "margin_scan",
            "case": case,
            "seed": 1,
            "margin_scan": {
                "error_model": {k: v for k, v in em.items() if k != "n_fractions"},
                "margins_mm": [2, 3, 4, 5, 6, 7, 8, 9],
                "constraint": {
                    "structure": "CTV",
                    "metric": "Dmin",
                    "op": ">=",
                    "value": 95.0,
                    "pct_of_rx": True,
                },
                "n_courses": 600,
            },
        },
    }


__all__ = ["ErrorModel", "demo_configs", "run_experiment"]
