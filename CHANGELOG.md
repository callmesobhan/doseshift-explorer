# Changelog

## 0.1.0 (unreleased)

First public release.

- Perturbation framework: physical (rigid shift/rotation, smooth deformation) and numerical
  (dose-grid offset, coarsening, interpolation, structure sampling), composable and pure.
- Exact weighted DVH, metric grammar (D95%, V20Gy, HI, CI, GI, gEUD), constraints.
- Monte Carlo of fractionated courses (sampled and closed-form `blur` engines); van Herk margin
  recipe check.
- Grid-assumption audit with setup-equivalent shift; measured `contour_variation`.
- DICOM-RT reader (CT/MR, RTDOSE, RTSTRUCT, RTPLAN prescription), header-reference case discovery,
  synthetic DICOM writer; SimpleITK registration.
- CLI, YAML experiments, Streamlit explorer, public-data downloader (`scripts/download_tcia.py`).
- Within-plan `resolution_scan` (does the cost of a half-voxel error scale with voxel size?),
  mutual-information CT-CBCT registration with ground-truth validation and a diagnostics script
  (`scripts/cbct_registration.py`).
- Three independent reviews acted on (engineering/release, scientific, registration/discovery/plots):
  loud registration failures, unambiguous case discovery, sign-aware equivalent shifts, consistent
  conformity index, hardened readers and downloader, corrected claims (see `docs/validation.md`).
