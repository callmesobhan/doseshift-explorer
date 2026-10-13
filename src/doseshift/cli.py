"""Command line interface: ``doseshift {demo, run, make-phantom, inspect, discover}``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

from . import __version__


def _cmd_demo(args) -> int:
    from .runner import demo_configs, run_experiment

    out = Path(args.out)
    for name, cfg in demo_configs(args.courses).items():
        print(f"[demo] {name} ...", flush=True)
        files = run_experiment(cfg, out / name)
        print("       " + ", ".join(sorted(p.name for p in files.values())))
    print(f"[demo] results in {out}/")
    return 0


def _cmd_run(args) -> int:
    from .config import load_config
    from .runner import run_experiment

    cfg = load_config(args.config)
    out = Path(args.out) if args.out else Path("outputs") / Path(args.config).stem
    files = run_experiment(cfg, out)
    for name, path in files.items():
        print(f"{name:20s} {path}")
    return 0


def _cmd_make_phantom(args) -> int:
    from .dicomrt import write_case
    from .synthetic import make_case

    paths = write_case(
        make_case(args.seed, args.ctv_margin), args.out, reverse_dose_z=args.reverse_dose_z
    )
    for k, v in paths.items():
        print(f"{k:9s} {v}")
    return 0


def _cmd_inspect(args) -> int:
    from .dicomrt import alignment_report

    print("\n".join(alignment_report(args.ct_dir, args.dose, args.rtstruct)))
    return 0


def _cmd_discover(args) -> int:
    from .dicomrt import case_config, discover_with_report
    from .state import ROLES

    if args.roles and not args.write_config:
        raise ValueError("--roles needs --write-config FILE (where to write the config)")
    roles = {}
    for item in (args.roles or "").split(",") if args.roles else []:
        name, sep, role = item.partition("=")
        if not sep or role not in ROLES:
            raise ValueError(f"bad --roles entry {item!r}: use ROI_NAME=ROLE with ROLE in {ROLES}")
        roles[name] = role
    cases, skipped = discover_with_report(args.patient_dir)
    for reason in skipped:
        print(f"skipped RTSTRUCT: {reason}")
    if not cases:
        print("no analysable case found (see the reasons above, if any)")
        return 1
    for i, c in enumerate(cases):
        print(f"case {i}: " + c.summary())
    if roles:
        if not 0 <= args.case < len(cases):
            raise ValueError(f"--case {args.case} out of range (found {len(cases)} case(s))")
        case = cases[args.case]
        cfg = case_config(case, roles, args.experiment, args.prescription)
        Path(args.write_config).write_text(yaml.safe_dump(cfg, sort_keys=False))
        print(f"wrote {args.write_config}")
    return 0


def _positive_int(text: str) -> int:
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return value


DISCLAIMER = "Research and teaching tool. Not validated for clinical use."


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="doseshift", description=__doc__, epilog=DISCLAIMER)
    p.add_argument("--version", action="version", version=f"doseshift {__version__}")
    p.add_argument("--debug", action="store_true", help="show full tracebacks on errors")
    sub = p.add_subparsers(dest="command", required=True)

    d = sub.add_parser("demo", help="run every experiment on the synthetic phantom")
    d.add_argument("--out", default="outputs/demo")
    d.add_argument("--courses", type=_positive_int, default=300, help="Monte Carlo courses")
    d.set_defaults(fn=_cmd_demo)

    r = sub.add_parser("run", help="run one experiment from a YAML config")
    r.add_argument("config")
    r.add_argument("--out", help="default: outputs/<config name>")
    r.set_defaults(fn=_cmd_run)

    m = sub.add_parser("make-phantom", help="write the synthetic phantom as DICOM-RT")
    m.add_argument("out")
    m.add_argument("--seed", type=int, default=0)
    m.add_argument("--ctv-margin", type=float, default=5.0)
    m.add_argument(
        "--reverse-dose-z",
        action="store_true",
        help="store RTDOSE frames in descending z (tests reader robustness)",
    )
    m.set_defaults(fn=_cmd_make_phantom)

    i = sub.add_parser("inspect", help="report dose-to-CT geometric alignment of a DICOM case")
    i.add_argument("--ct-dir", required=True)
    i.add_argument("--dose", required=True)
    i.add_argument("--rtstruct")
    i.set_defaults(fn=_cmd_inspect)

    c = sub.add_parser("discover", help="pair image/RTSTRUCT/RTPLAN/RTDOSE in a downloaded folder")
    c.add_argument("patient_dir")
    c.add_argument("--roles", help="e.g. TV=ptv,Cochlea=oar; with --write-config makes a config")
    c.add_argument("--case", type=int, default=0, help="which discovered case (default 0)")
    c.add_argument(
        "--experiment", default="grid_audit", choices=["grid_audit", "sweep", "resolution_scan"]
    )
    c.add_argument("--prescription", type=float, help="Gy; needed when the data has no RTPLAN")
    c.add_argument("--write-config", metavar="FILE", help="write a ready-to-run YAML here")
    c.set_defaults(fn=_cmd_discover)

    args = p.parse_args(argv)
    try:
        return args.fn(args)
    except (OSError, ValueError, KeyError, RuntimeError, yaml.YAMLError) as exc:
        if args.debug:
            raise
        print(f"error: {exc}", file=sys.stderr)
        print("(re-run with --debug for the full traceback)", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
