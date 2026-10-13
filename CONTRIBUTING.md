# Contributing

Thanks for looking. This is a research tool; correctness and reproducibility matter more than
features.

- **Never commit patient data.** No DICOM files, no image slices, no screenshots of real anatomy.
  `data/*` (except its README) and `outputs/` are git-ignored; do not force-add. Plots of numbers
  derived from public data are fine, with attribution (see `data/README.md`).
- Set up: `pip install -e ".[all]"`, then `pytest` and `ruff check . && ruff format --check .`.
  Ruff is pinned (see the `dev` extra) so formatting does not drift.
- A new perturbation or metric needs: a test that its identity setting changes nothing, and a line
  in `docs/conventions.md` if it defines a convention. See *Extending* in `docs/design.md`.
- Claims in the docs should be reproducible from the code. If you add a number to a README or
  `docs/experiments.md`, say which command produces it.
- Real-data tests (`tests/test_real_data.py`) skip when `data/raw/` is absent; fetch the data with
  `scripts/download_tcia.py` (dry run by default).

This software is for research and teaching only and is **not validated for clinical use**.
