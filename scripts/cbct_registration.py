#!/usr/bin/env python3
"""CT <-> CBCT rigid registration diagnostics for a Pancreatic-CT-CBCT-SEG patient.

    python scripts/cbct_registration.py \
        data/raw/pancreas-ct-cbct-seg/Pancreas-CT-CB_014 --out cbct_014.csv

What it does, per CBCT series
1. Pairs each original CBCT with the dataset's "Aligned CT" series (same content placed into the
   planning frame) by same-modality registration error, and reads off the dataset's own placement.
2. Registers the CBCT to the planning CT by mutual information from several starting translations
   (the dataset placement plus offsets), and reports how much the answer depends on the start.

What it can and cannot tell you
* The CBCT frame is unrelated to the planning frame and the collection has no RTPLAN (no isocentre),
  so a translation from this registration is **not** a measured setup error. The informative
  numbers are the rotation, the spread over starts (is the optimum well defined?), and the
  difference from the dataset's placement.
* There is no ground truth for the real pairs. The mutual-information path is validated separately
  on known motions across scanner-like differences (tests/test_registration.py and
  tests/test_real_data.py); a small start spread on real data supports, but does not prove, a
  result. Start spreads of several mm mean the optimum is not well defined: do not use that pair.
* Research use only; not validated for clinical image guidance.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from doseshift.dicomrt import discover_cases, load_image_series
from doseshift.registration import intensity_centroid_mm, register_rigid

START_OFFSETS_MM = [
    (0, 0, 0),
    (12, 0, 0),
    (-12, 0, 0),
    (0, 12, 0),
    (0, -12, 0),
    (0, 0, 12),
    (0, 0, -12),
]


def find_series(patient: Path):
    cases = discover_cases(patient)
    planning = next((c for c in cases if any(r.endswith("planCT") for r in c.roi_names)), None)
    if planning is None:
        sys.exit("no planning-CT structure set (ROI ending in 'planCT') found")
    aligned_dirs = {c.image_dir for c in cases if c is not planning}
    plan_grid, plan_vals, plan_for = load_image_series(planning.image_dir)
    cbct = {}
    for d in sorted(p for p in patient.iterdir() if p.is_dir() and p.name.startswith("CT_")):
        if d == planning.image_dir or d in aligned_dirs:
            continue
        g, v, f = load_image_series(d)
        if f != plan_for:  # own frame of reference: an original CBCT
            cbct[d] = (g, v)
    aligned = {d: load_image_series(d)[:2] for d in aligned_dirs}
    return (plan_grid, plan_vals), aligned, cbct


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("patient_dir", type=Path)
    ap.add_argument("--out", type=Path, help="write the per-start results as CSV")
    args = ap.parse_args(argv)
    planning, aligned, cbct = find_series(args.patient_dir)
    pg, pv = planning

    # 1. the dataset's placement: same content, so mean squares is appropriate
    placement = {}
    for ad, (ag, av) in aligned.items():
        scored = []
        for cd, (cg, cv) in cbct.items():
            t0 = intensity_centroid_mm(cg, cv) - intensity_centroid_mm(ag, av)
            r = register_rigid(
                (ag, av),
                (cg, cv),
                metric="mse",
                initial_translation_mm=tuple(t0),
                shrink_factors=(4, 2),
                smoothing_sigmas_vox=(2.0, 1.0),
            )
            scored.append((r.metric_value, cd, r))
        err, cd, r = min(scored, key=lambda x: x[0])
        placement[cd] = (ad.name, np.array(r.motion.translation_mm), r.motion.rotation_deg, err)

    rows = []
    for cd, (cg, cv) in cbct.items():
        ad_name, t_ds, rot_ds, mse = placement[cd]
        for k, off in enumerate(START_OFFSETS_MM):
            r = register_rigid(
                (pg, pv),
                (cg, cv),
                metric="mi",
                initial_translation_mm=tuple(t_ds + np.array(off, float)),
                seed=k,
                shrink_factors=(4, 2),
                smoothing_sigmas_vox=(2.0, 1.0),
            )
            t = np.array(r.motion.translation_mm)
            rows.append(
                {
                    "cbct": cd.name,
                    "aligned_ct": ad_name,
                    "pair_mse": mse,
                    "start": k,
                    "start_offset_mm": str(off),
                    **{f"t_{a}": float(v) for a, v in zip("xyz", t, strict=True)},
                    **{
                        f"dt_vs_dataset_{a}": float(v) for a, v in zip("xyz", t - t_ds, strict=True)
                    },
                    **{
                        f"rot_{a}_deg": float(v)
                        for a, v in zip("xyz", r.motion.rotation_deg, strict=True)
                    },
                    "mi_metric": r.metric_value,
                    "stop": r.stop_condition,
                }
            )
    df = pd.DataFrame(rows)
    if args.out:
        df.to_csv(args.out, index=False)

    print("Dataset placement of each CBCT into the planning frame (same-modality registration):")
    for cd, (ad_name, t_ds, rot_ds, mse) in placement.items():
        print(
            f"  {cd.name} -> {ad_name}: translation {np.round(t_ds, 2)} mm, rotation "
            f"{np.round(rot_ds, 2)} deg, residual MSE {mse:.0f}"
        )
    print("\nMutual-information registration to the planning CT, over start offsets:")
    for cd in cbct:
        d = df[df["cbct"] == cd.name]
        dt = d[[f"dt_vs_dataset_{a}" for a in "xyz"]].to_numpy()
        spread = d[[f"t_{a}" for a in "xyz"]].std().max()
        rot = d[[f"rot_{a}_deg" for a in "xyz"]].median().to_numpy()
        med = np.median(dt, axis=0)
        print(
            f"  {cd.name}: median difference from dataset placement {np.round(med, 1)} mm "
            f"(|d| = {np.linalg.norm(med):.1f} mm); start spread (max axis std) {spread:.2f} mm; "
            f"median rotation {np.round(rot, 2)} deg"
        )
    print(
        "\nReminder: translations here are NOT setup errors (unrelated frames, no isocentre in\n"
        "the data). A small start spread means a well-defined optimum, not a correct one."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
