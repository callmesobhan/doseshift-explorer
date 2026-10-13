# Validation

What is checked, against what, and where. `pytest` runs all of it (~40 s; 183 tests, 5 of which need the downloaded real data and skip without it).

| Claim | Checked against | Test |
|---|---|---|
| Moving the patient by +t equals locating the dose at −t | itself, two code paths (the core invariant) | `test_perturb.py::test_patient_shift_equals_opposite_dose_grid_shift` |
| Identity perturbations change nothing; input state is never mutated | exact equality | `test_zero_perturbations_are_identity`, `test_state_is_not_mutated` |
| Rotation convention is right-handed `Rz Ry Rx`; inverse undoes forward | closed-form cases | `test_geometry.py` |
| Motion sign, rotation order and centre handling match an **independent** SimpleITK resampler | moving image made by SimpleITK (not the package's own warp), registered back, on oblique, negative-k, left-handed and oblique+negative-k grids with an off-centre rotation: point error < 0.3 mm when the whole body is in the field of view (reviewer: 0.002-0.084 mm); 0.5-1.2 mm when the field truncates the body | `test_registration.py` |
| Cross-scanner mutual information (different grid, frame origin, intensity mapping) | known motion recovered to ≤ 0.04 mm / 0.5° (phantom) and ≤ 1.0 mm / 0.14° (real pancreas CT) from starts up to 12 mm off | `test_registration.py`, `test_real_data.py` |
| DVH metrics are exact (no binning) | hand-computed weighted cases, monotonicity, extremes | `test_dvh_metrics.py` |
| Paddick CI | hand-computed cases; independence from dose-grid spacing on a plan whose isodose is the target surface (1 mm: 1.00, 5 mm: > 0.88) | `test_paddick_ci_perfect_and_partial`, `test_conformity_index_is_nearly_independent_of_dose_grid_spacing` |
| `GI`/`PIV` are invariant to setup shift but not to grid assumptions | direct comparison | `test_dose_cloud_metrics_ignore_setup_shift_but_not_grid_assumptions` |
| DICOM reader geometry (including **reversed** `GridFrameOffsetVector`) | write → read round trip: CT within 0.5 HU, dose interpolated identically, RTSTRUCT masks identical | `test_dicomrt.py` |
| `blur` engine is the many-fraction limit of `sampled` | 1500-fraction run vs closed form | `test_blur_engine_matches_sampler_in_many_fraction_limit` |
| Monte Carlo is **consistent with the van Herk margin recipe** under the phantom's assumptions | P(CTV Dmin ≥ 95 % Rx) = 0.89 at 6.95 mm (3D) and ≈ 0.9 at the 1D recipe margin; the test accepts 0.84-0.95 | `test_reproduces_van_herk_margin_recipe` |
| Equivalent shift is monotone, sign-aware and a lower bound; a 3.0 mm (one-voxel) dose offset is recovered at the right sign for rectum D2% and never exceeded for PTV D95 | synthetic offsets with known displacement | `test_experiments.py`, `test_review_fixes.py` |
| Example configs stay valid; every experiment writes tables, figures, provenance; DICOM → CLI → experiment works | end-to-end | `test_runner_cli.py` |
| The Streamlit app runs and reacts to a shift | `streamlit.testing.AppTest` | `test_app.py` |

## Real data

Run against public clinical data (`tests/test_real_data.py`, skipped when `data/raw/` is absent):

| Check | Result |
|---|---|
| Reader on RTDOSE with a **negative** `GridFrameOffsetVector` (VS-SEG, dz = −1 mm) | loads; k-step is negative in the affine |
| Reader on dose grids unrelated to the CT grid (pancreas: 2.5 × 2.5 × 3 mm and 1.25 × 1.25 × 1.5 mm vs 1.17 × 1.17 × 3 mm CT) | loads; sub-voxel phase reported by `inspect` |
| Plan sanity: tumour D95% vs the RTPLAN's 13 Gy prescription | 13.2-13.3 Gy |
| Grid-audit scale vs voxel size | ½-voxel offset ≈ 0.5-0.8 mm on 1 mm grids, 0.5-2 mm on 2.5 mm grids |
| Pairing: `discover` vs hand-picked files | found that hand-picking had chosen the wrong RTSTRUCT and RTDOSE for a patient |

Real data exposed three defects in the first version, all fixed: the reader accepted only CT (the
Gamma Knife set is MR); nothing checked that an RTSTRUCT referenced the loaded image series; and a
skipped-contour warning did not say which structure. The reader had previously been tested only on
DICOM written by this repo's own writer.

## What validation does not cover

- Agreement with a commercial TPS or with measured data. The phantom's dose is analytic.
- Agreement of results with an independent dosimetry or with measured setup errors; real-data
  tests check plausibility and internal consistency, not clinical truth.
- Anything outside the static-dose-cloud regime ([limitations](limitations.md)).

## Independent review before publication

Two independent reviewers (separate agents given the code and no conclusions) were asked to find
problems; a first attempt at each was lost to an API rate limit and the scientific review was
re-run on a private snapshot of the code. Their reports are model output and were treated as leads:
every finding below was re-checked before being acted on.

**Engineering and release review** (fixes covered by `tests/test_hardening.py`,
`tests/test_download_script.py`):

| Finding | Fix |
|---|---|
| Dataset attribution wrong (credited the methods paper's authors) and incomplete | exact TCIA citations from DataCite, CC BY link, "changes made" statement, `CITATION.cff` |
| Downloader treated a corrupt or truncated download as complete and left a zip behind | atomic extract into `*.partial`, `.complete` marker only after zip check and file-count match, retries, manifest after success, dry run creates nothing |
| Silent wrong numbers: duplicate ROI names, `DoseGridScaling` 0/NaN/negative, prescription 0/NaN, non-planar contours dropped | error / warning naming the file and the problem |
| Raw tracebacks or opaque errors on malformed DICOM (NaN contours, missing tags, bad ROI refs, mixed slice sizes, undecodable pixel data) | `ValueError` naming file and attribute; decoder hint |
| Structure partly or wholly outside the dose grid read 0 Gy silently; empty structures | warning with the percentage; error naming the structure |
| `discover` crashed on a non-DICOM `*.dcm` and missed extension-less files | content sniffing |
| CLI tracebacks, `--courses 0`, `n_courses: 0`, relative paths from the wrong directory | clean `error:` messages, `--debug` for tracebacks |
| Real-data runs rendered patient anatomy into `geometry.png` | off by default for DICOM cases |
| Provenance lacked matplotlib/pydicom/SimpleITK versions, platform, input identity; `provenance.json` could not be re-run | all added; a provenance file is a valid config |
| `GridFrameOffsetVector[0]` ignored | origin = IPP + normal x offset[0] |
| `equivalent_shift_mm` shown for organs whose volume changed too much; one noise floor for every unit | NaN for those rows; dimensionless indices use a smaller floor |
| Unpinned ruff could turn CI red; editable install only; no least-privilege CI | pinned ruff, non-editable test install, build job, `permissions: contents: read` |

Reproducibility, as measured by the reviewer: two `doseshift demo` runs on one machine give
byte-identical CSVs and PNGs; across Python/numpy/matplotlib versions CSVs agree to about 1e-13 and
the PNGs differ in pixels. Results are therefore **numerically** reproducible, not bit-for-bit across
environments. `provenance.json` records what is needed to rerun.

**Scientific and numerical review** (fixes covered by `tests/test_review_fixes.py` and
`tests/test_hardening.py`). The review also confirmed things worth knowing: grid index/world maps, the
DICOM affine (against an independent PS3.3 construction and SimpleITK, with 0.5 x 1.25 mm anisotropic
pixel spacing, 20° obliquity and both frame-offset directions), `RigidMotion`, and the DVH (5421
brute-force comparisons, 0 mismatches) are correct, and every number in the real-data tables
reproduced.

| Finding | Fix or change of claim |
|---|---|
| "MC reproduces van Herk (0.89), not a coincidence" is not robust: 0.79-0.91 across one-at-a-time setup changes; `blur` ignores `n_fractions` | claim reworded to a consistency check; sensitivity documented (experiments R4); binomial standard errors reported; margin scan marked phantom-only |
| `equivalent_shift_mm` matched magnitude only, so a 3.0 mm offset was reported as 1.6 mm; for offsets the "equivalent" is just the displacement; "scales with voxel size" was true by construction | sign-aware matching; actual displacement reported beside it; documented as a lower bound; R5 rewritten; scaling claim withdrawn |
| `StructureSampling` silently dropped the part of a structure outside the dose grid (a +58 Gy fake "numerical" change) | out-of-grid voxels are kept (0 Gy, as native); volume conserved; warning when a structure is outside the dose grid |
| Paddick CI mixed an interpolated TV_PIV with a node-count PIV, biasing CI with grid spacing | both from the same interpolated field on a refined lattice |
| `blur` under-blurred when sigma < one dose voxel (7× at 1 mm on a 3 mm grid) and applied sigma per index axis | grid refined first; sigma mapped onto world axes; oblique grids rejected |
| `_dice` across different grids was asymmetric and could exceed 1 | symmetric estimator, capped at 1 |
| Cochlea row was a quantisation artefact (5 samples left) | warning below 30 samples; row withdrawn as a conclusion |
| `DoseGridOffset(v)` equals a shift of `-A v`, not `-v * spacing`; audit labels were world axes but are index axes | docs corrected; labels are i, j, k |
| Smaller edge cases: `D{x}cc` on a smaller structure returned Dmin, NaN dose gave finite D95, `gEUD(a=0)` crashed, `SmoothDeformation` realised 0.65-1.25 mm for a request of 2, candidate box too small for coarse structure voxels, coarsening dropped trailing voxels silently | NaN / error / geometric mean, zero-mean deformation with the requested realised RMS, margin from voxel extent, warning |

Not covered by that review: the downloader (engineering review) and the plots, app, `registration.py`,
`discover.py` and the timing figures (third review, below).

**Third review: registration, case discovery, plots, app, timing claims** (fixes covered by
`tests/test_registration.py`, `tests/test_discover_adversarial.py`, `tests/test_review_fixes.py`,
`tests/test_real_data.py`). It confirmed that the registration geometry, the plot orientation labels
("left on the right of the image", "anterior on the left", checked against the LPS phantom), the
signed sweep concatenation and the Monte Carlo band direction are correct.

| Finding | Fix or change of claim |
|---|---|
| `register_rigid` returned an **identity motion with no error** when the images did not overlap (also on the real pancreas CT vs CBCT, whose frames are unrelated): "patient exactly on plan" | `ValueError` when the foreground barely overlaps at the start (with the fix to try), `RuntimeError` if the metric is never evaluated, NaN and thin-volume (< 4 voxels) guards; result carries `converged`, `overlap_fraction` and per-level iterations; warns on the iteration limit |
| Discovery silently took the first of several plans, doses or image series, picked beam doses over plan doses depending on file names, took the max of several target doses, and dropped unplaceable RTSTRUCTs | one case per plan and per image series; plan-level physical dose preferred, ties not paired; prescription withheld for several targets or phases; images grouped by series UID; skipped RTSTRUCTs reported with reasons (CLI explains why nothing was found); dose frame-of-reference checked |
| App printed "nan mm" and said "magnitude-matched" while computing a signed match | NaN-aware wording ("no rigid shift can reproduce this"), "lower bound", `< 0.1 mm` rule |
| More than 8 structures aborted the whole run | neutral grey beyond 8 with a warning |
| Colours were not stable per structure across figures; sweeps reused structure colours; margin and histogram plots were blue regardless of structure | one colour order per plan passed to every figure; sweep axes use their own slots; bars, histogram and margin scan take the structure's colour |
| Subtitles clipped; "5-95 % of courses" band mislabelled; threshold had no unit; sweep legend listed the positive end first; `use_container_width` deprecated | shortened, "pointwise percentile", units, "+ = patient left" wording, argument dropped |
| Registration accuracy quoted as a general figure, though it holds for full-field images only | scope stated here and in `docs/experiments.md` R8; truncated-field test added |
| Timing claims wrong (3 structures / 15 k points) | re-measured: 4 structures, ~17 k points; `docs/design.md` corrected |

## Findings from building it

Two checks were corrected during development and are worth keeping in mind:

- *D98 vs minimum dose.* A first margin-scan used CTV D98, and the plan appeared to need much less
  margin than the recipe predicts (1D error, 3 mm margin: P = 0.85 with D98 vs 0.52 with `Dmin`).
  That is correct behaviour, not a bug: for an ellipsoid, D98 discards a cap a few mm deep, so the
  criterion is weaker than the recipe's minimum-dose criterion. With `Dmin` the recipe is
  reproduced.
- *Recipe factor.* 2.5 in the recipe is the 3D value. A 1D error model needs 1.28 (exposed as
  `ErrorModel.van_herk_margin_mm(systematic_factor=...)`).
