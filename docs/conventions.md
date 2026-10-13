# Conventions

Geometry conventions are where this kind of code silently goes wrong, so they are fixed in one
place and each is covered by a test.

## Coordinates and arrays

| | Convention |
|---|---|
| World frame | DICOM patient coordinates (LPS), mm. **+x** patient left, **+y** posterior, **+z** superior |
| Vector order | `(x, y, z)` everywhere (points, spacing, shifts, rotations) |
| Array order | `[k, j, i]` = `[z, y, x]` = (slice, row, column), as pydicom / SimpleITK return them |
| `Grid.affine` | 4×4, continuous index `(i, j, k, 1)` → world `(x, y, z, 1)` |
| Voxel position | `affine @ (i, j, k)` is the **centre** of voxel `(i, j, k)` |
| Outside the dose grid | dose = 0 Gy (no extrapolation) |

## DICOM → `Grid`

| DICOM | Used as |
|---|---|
| `ImagePositionPatient` | centre of the first voxel (affine translation) |
| `ImageOrientationPatient` | `[row_dir, col_dir]`; `i` runs along `row_dir`, `j` along `col_dir` |
| `PixelSpacing` | `[row spacing, column spacing]`: `i`-step = `row_dir × PixelSpacing[1]`, `j`-step = `col_dir × PixelSpacing[0]` |
| CT slice position | sorted by projection on `row_dir × col_dir`; spacing must be uniform |
| RTDOSE `GridFrameOffsetVector` | slice offsets from `ImagePositionPatient` along the normal (first entry included). **A decreasing vector is honoured**: the `k` step is negative in the affine (tested with `write_case(..., reverse_dose_z=True)`) |
| RTDOSE values | `pixel_array × DoseGridScaling`; only `DoseUnits = GY` accepted |
| RTSTRUCT contours | planar polygons, one contour → one slice; a voxel is inside if its **centre** is inside the polygon; contours on the same slice combine by XOR, so nested contours make holes |

Image series may be CT (HU, rescale applied) or MR (raw intensities). An RTSTRUCT's contours belong
to the image series named in its `ReferencedFrameOfReferenceSequence`; `inspect` and
`load_plan_state` check this, and `discover` uses it to pair files. RTDOSE belongs to the RTPLAN it
references, which references the RTSTRUCT it was planned on.

`doseshift inspect` prints the fractional position of the dose grid's first voxel in CT index space
(the sub-voxel phase between the two grids), FrameOfReferenceUID agreement, and extents.

## Patient motion

`RigidMotion(translation_mm, rotation_deg, center_mm)` describes where the *patient's anatomy*
ends up: a planned point `p` is found at

```
p' = R (p − c) + c + t        R = Rz(rz) · Ry(ry) · Rx(rx)    (right-handed, about fixed world axes)
```

and the dose is evaluated at `p'`. The dose stays fixed in the room. Rotation centre defaults to
the PTV centroid. In the audit, `DoseGridOffset(v)` *locates the dose* `v` voxels further along its own
**index** axes (i, j, k), which displaces it by `A·v` in the room (`A` = the 3×3 part of the grid
affine) and is equivalent to patient motion `t = −A·v`. Only for an axis-aligned grid with positive
steps is that `−v · spacing` along x, y, z. For the Gamma Knife doses (slice direction negative) a
"+½ voxel (k)" offset moves the dose towards −z. This is why audit scenarios are labelled by index
axis (i, j, k), not by world axis. (Verified against an oblique, negative-k grid by an independent
reviewer.)

### Crosswalk with SimpleITK / ITK

Registering the treatment image (moving) to the planning CT (fixed) returns `T` with
`T(p_fixed) = p_moving`, i.e. exactly `RigidMotion`. `registration.register_rigid` reads the
transform's 3×3 matrix, translation and centre and re-expresses the rotation as `Rz Ry Rx` Euler
angles, so ITK's own Euler parameterisation does not matter. `tests/test_registration.py` checks the
sign and the rotation order by injecting a known motion and recovering it.

## Metrics

All metrics are computed on weighted samples (`PointCloud.volumes`, mm³), not on a histogram.

| Spec | Definition | Unit |
|---|---|---|
| `Dmean`, `Dmax`, `Dmin` | weighted mean / max / min sampled dose | Gy |
| `D{x}%` | highest dose D such that ≥ x % of the volume receives ≥ D | Gy |
| `D{x}cc` | same, for x cc of absolute volume | Gy |
| `V{x}%` | % of volume receiving ≥ x % of the prescription | % |
| `V{x}Gy` | % of volume receiving ≥ x Gy | % |
| `V{x}Gy_cc` | volume receiving ≥ x Gy | cc |
| `gEUD(a=n)` | (Σ vᵢ dᵢᵃ)^(1/a); a ≫ 1 serial organs, a < 0 tumours | Gy |
| `HI` | ICRU 83: (D2% − D98%) / D50% | – |
| `CI` | Paddick: TV_PIV² / (TV × PIV) | – |
| `GI` | V50%Rx / V100%Rx, volumes of the interpolated dose field | – |

`D2%` is "near-maximum"; `Dmax` is the single hottest sample and is noisy by construction.
`D{x}cc` is NaN for a structure smaller than x cc; `gEUD(a=0)` is the geometric mean; non-finite dose is
rejected. Structures with fewer than 30 sample points trigger a warning (metrics are then quantisation).
`PIV` (volume where the dose is ≥ Rx) and the `GI` volumes use the whole dose grid and do not move
with the patient; `TV_PIV` is the part of the (moved) target that receives ≥ Rx. Both come from the
**same interpolated field**: PIV is evaluated on a lattice refined to about 1 mm (cached per dose grid), not by
counting dose-grid nodes, because mixing the two biased CI with the dose-grid spacing (CI of a plan whose
isodose is exactly the target surface was 0.97 / 0.91 / 0.82 / 0.58 on 1 / 2 / 3 / 5 mm grids; it is now
1.00 / 0.98 / 0.98 / 0.92). A grid too coarse to represent the penumbra still lowers CI; that is real.

References: HI follows ICRU Report 83 (*J. ICRU* 10(1), 2010); CI is Paddick, *J Neurosurg* 93
Suppl 3:219-222 (2000), [doi:10.3171/jns.2000.93.supplement_3.0219](https://doi.org/10.3171/jns.2000.93.supplement_3.0219);
the margin recipe is van Herk, Remeijer, Rasch & Lebesque, *Int J Radiat Oncol Biol Phys* 47(4):1121-1135
(2000), [doi:10.1016/S0360-3016(00)00518-6](https://doi.org/10.1016/S0360-3016(00)00518-6).

Keys in tables and constraints are `STRUCTURE:metric`, e.g. `PTV:D95%`.

## Setup-error model

`ErrorModel(sys_sigma_mm, rand_sigma_mm, sys_sigma_deg, rand_sigma_deg, n_fractions)`:
`Σ` (systematic) is drawn once per simulated course, `σ` (random) per fraction, axes (x, y, z).
Delivered dose = mean over fractions of the planned dose at the displaced points. Margin recipe:
`M = a·Σ + 1.64·(√(σ² + σp²) − σp)` with `σp = 3.2 mm`, `a = 2.5` (3D) or `1.28` (1D).

## Equivalent rigid shift

`experiments.equivalent_shift_mm(sweeps, key, delta)` inverts translation sweeps along ±x, ±y, ±z to find
the smallest shift that changes a metric by `delta`. It is **sign-aware** (a shift must move the metric the
same way; NaN if none can) and a **lower bound** on any displacement behind the change, because only
axis-aligned directions are searched and the metric may be steeper along one of them than along the true
displacement. For `DoseGridOffset` scenarios the true displacement `|A v|` is reported alongside
(`displacement_mm`). 0 means below the noise floor (`min_delta`, 0.1 in Gy / cc / %, one twentieth of
that for dimensionless indices), inf means beyond the swept range.

## Blur engine constraints

`montecarlo.simulate(method="blur")` convolves the planned dose with the random-error Gaussian. It
needs a dose grid whose axes are a signed permutation of the patient axes (it raises otherwise) and maps
sigma onto the grid's index axes accordingly; if sigma is below one dose voxel on an axis the dose is first
refined linearly so the filter is resolved (filtering at sub-voxel sigma on the original grid under-blurred
about 7× at σ = 1 mm on a 3 mm grid). It is the infinite-fraction limit and ignores `n_fractions`.
