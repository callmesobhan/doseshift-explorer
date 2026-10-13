# Roadmap

## Done (v0.1)

Core model and perturbations, exact DVH and metric grammar, Monte Carlo with two engines, grid
audit with setup-equivalent shift, DICOM-RT reader/writer, SimpleITK registration, CLI, experiments,
Streamlit explorer, 65 tests (144 after the real-data work and the independent reviews below).

## Done since v0.1

Real data: two public TCIA collections (5 patients), `scripts/download_tcia.py`, header-reference
case discovery, MR support, `contour_variation`, and the real-data results in
[experiments.md](experiments.md).

Independent engineering review before publication (hostile-input handling, downloader integrity,
attribution, provenance, CI, ignore rules), see [validation.md](validation.md).

## Next (v0.2)

- **Measure setup error from the pancreas CBCTs.** The original CBCT series (own frame of
  reference, 0.51 mm in-plane) are already downloaded; register them to the planning CT with
  `registration.register_rigid` (mutual information, since CT and CBCT intensities differ) and fit
  an empirical `ErrorModel` for R3.
- **Equivalent shift along a known direction.** Report the equivalent along the displacement
  direction (and a min/median/max over many directions), not only the best axis-aligned one.
- **Pooled analysis across patients** so that voxel-size scaling can actually be tested.
- **More patients / a head-neck set.** Add pancreas and VS-SEG patients as disk allows; request
  HMM-RT (gated) or the TCIA collections that need an agreement, if you want them.

- **Contour-level rasterisation perturbation.** Re-rasterise RTSTRUCT contours at sub-voxel offsets
  and with partial-volume rules, so the one convention not yet perturbed (mask sampling from
  contours) joins the audit. Needs the contours kept on `MaskStructure`.
- **Gamma / dose-difference between perturbed dose grids** as a second lens next to DVH metrics.
- **Real-case notebook** on a public collection (see [data/README.md](../data/README.md)): inspect,
  audit, compare to the phantom's ranking.
- **Empirical error model from registration:** batch `register_rigid` over repeat scans and fit Σ, σ.

## Later

- Deformable registration (SimpleITK B-spline/demons) to replace `SmoothDeformation` with measured
  fields; deformable dose accumulation.
- Radiobiology: EQD2 conversion, TCP/NTCP wrappers, gEUD-based constraints in the Monte Carlo.
- Pluggable dose engines (re-compute instead of static cloud) for heterogeneity / proton cases.
- RTPLAN reader (beam geometry, isocentre from the plan instead of the PTV centroid).
- Robust-optimisation comparison: nominal vs margin vs probabilistic plans under the same error model.
