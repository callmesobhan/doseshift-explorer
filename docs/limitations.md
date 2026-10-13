# Limitations

Read before using any number from this repo for anything beyond research and teaching.

**Physics**

- **Static dose cloud.** Dose is not recomputed after the patient moves. Valid to first order for
  photons in near-homogeneous tissue; invalid for protons or where shifts change the radiological
  path (lung, bone, couch, air cavities) and for IMRT/VMAT interplay effects.
- **Rigid + synthetic smooth deformation only.** `SmoothDeformation` is a controlled stand-in, not
  a model of real anatomical change.
- **Dose is evaluated per structure point and fractions are averaged**, with no radiobiological
  weighting (no EQD2, no TCP/NTCP).

**Numerics**

- **The equivalent rigid shift is a lower bound**, sign-aware, over six axis-aligned directions. It is
  not an estimate of an error and not a measurement of setup error; for dose-grid offsets the true
  displacement is known and reported.
- **Paddick CI is a property of the dose representation.** CI and GI use the interpolated dose field;
  a grid too coarse to represent the penumbra still lowers CI (5 mm grid, 3 mm penumbra: about 0.92).
- **`blur` Monte Carlo** is the infinite-fraction limit, needs a dose grid aligned with the patient axes,
  and refines the dose when sigma is below one voxel; it ignores `n_fractions`. Coverage probabilities
  carry binomial standard errors; at the recipe margin the margin-scan probability ranges over
  0.79-0.91 as mask resolution, dose spacing or plateau level change.
- Dose beyond the outermost voxel *centres* is 0 Gy (the dose starts at the centre of the edge voxel,
  not its face); coarsening drops trailing voxels and shifts that edge inward.
- A structure with very few sample points (< 30) has quantised metrics; the tool warns.

- Dose outside the dose grid is 0 Gy.
- Structures are binary at the CT resolution unless `StructureSampling` is used; contour
  rasterisation uses voxel-centre-in-polygon, which is itself a convention that is not yet
  perturbed (see roadmap).
- DVH of a coarse structure is only as good as its point sample: `Dmin`, `Dmax` and small-volume
  `D0.03cc` are voxel-quantised and noisy.
- Cubic interpolation can overshoot near steep gradients (it is included to show exactly that).
- The setup-equivalent shift is magnitude-matched and direction-agnostic; it is a communication
  device, not an inverse problem solution.

**Real data**

- Two collections, five patients, two of each analysed in depth: the real-data results are
  illustrations of the method, not population statistics.
- The pancreas data has no RTPLAN and no PTV; the reference dose (ROI median) is derived from the
  dose itself and is not a prescription. Rx-dependent metrics are omitted there.
- `contour_variation` conflates organ motion, filling, residual setup error and *differences in
  contouring convention*; `rigid_comparable` only catches the volume part of that. A rigid
  equivalent is a communication device, not a measurement of setup error.
- Contour overlap here compares planning vs CBCT-day structure sets in a shared frame of reference;
  how well that frame was aligned is not verified.
- MR intensities are not calibrated (no HU); Gamma Knife frame/coordinate details are taken from
  the standard DICOM geometry only.

**Registration**

- Rigid registration is accurate to < 0.3 mm when the whole body is in the field of view and
  0.5-1.2 mm when the field truncates it (the case for CBCT). The returned translation is the
  displacement of the chosen rotation centre; pass the isocentre if you want a shift at the isocentre.
- Cross-scanner (mutual-information) registration has no ground truth on real pairs. Accept a result
  only if the optimum is start-independent (`scripts/cbct_registration.py` reports the spread); on the
  pancreas CBCTs that held for 1 of 4 pairs, and the frames are unrelated and carry no isocentre, so
  none gives a setup error. A registration that hits the iteration limit warns and sets
  `converged=False`; one that cannot start raises.

**Case discovery**

- `discover` pairs files from DICOM references and reports what it cannot place, but it cannot know
  your intent: with several plans or image series it returns one case each, and with a prescription
  that differs per target or phase it returns none (supply `prescription_gy`).

**Data and software**

- The phantom is analytic and prostate-like; effect sizes depend on dose gradient, grid spacing
  and structure geometry, so they will differ on real plans.
- The synthetic DICOM writer is minimal and not conformant enough to import into a clinical system.
- The DICOM reader assumes uniform slice spacing, absolute-Gy RTDOSE and planar contours; it reads
  CT and MR. RTPLAN is read only for the prescription and fraction count (not beam geometry).
- Not a medical device; not validated for clinical decision-making.
