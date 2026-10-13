# Data

**Never commit patient data.** `data/raw/`, `data/phantom/`, `outputs/` and `*.dcm` are git-ignored.
Respect each collection's licence and citation terms.

## Real data used here (downloaded, public, CC BY 4.0)

Both collections come from The Cancer Imaging Archive (TCIA); licence and file sizes below were read
from TCIA's own API. They are the only two collections in TCIA's *public* listing that carry
RTDOSE together with RTSTRUCT (checked on 2026-10-13 by scanning the modalities of all 156 public
collections; the listing changes over time).

| | Pancreatic-CT-CBCT-SEG | Vestibular-Schwannoma-SEG |
|---|---|---|
| Objects | CT (planning, 2 CBCT, 2 aligned), RTSTRUCT ×3, RTDOSE | MR (T2 / contrast T1), RTSTRUCT ×2, RTDOSE ×2, RTPLAN ×2 |
| Treatment | ablative pancreatic RT, MSKCC, deep-inspiration breath-hold | single-fraction Gamma Knife radiosurgery, 13 Gy |
| Size | 40 patients, 14.3 GB (231-595 MB / patient) | 242 patients, 27.7 GB (36-145 MB / patient) |
| Dose grid | 2.5 × 2.5 × 3 mm and 1.25 × 1.25 × 1.5 mm (patients 014, 017) | 1 mm isotropic, **negative** `GridFrameOffsetVector` |
| Data licence | CC BY 4.0 | CC BY 4.0 |
| Dataset (authors) | Hong et al. 2021, v2, [doi:10.7937/TCIA.ESHQ-4D90](https://doi.org/10.7937/TCIA.ESHQ-4D90) | Shapey et al. 2021, v2, [doi:10.7937/TCIA.9YTJ-5Q73](https://doi.org/10.7937/TCIA.9YTJ-5Q73) |
| Related paper | Han et al., *Med Phys* 2021 (methods paper; different author list from the dataset), [doi:10.1002/mp.14906](https://doi.org/10.1002/mp.14906): free author manuscript on PubMed Central, not a CC-licensed open-access article | Shapey et al., *Sci Data* 2021, [doi:10.1038/s41597-021-01064-w](https://doi.org/10.1038/s41597-021-01064-w): **open access, CC BY** |

**Downloaded for this repo (576 MB):** `Pancreas-CT-CB_014`, `Pancreas-CT-CB_017`, `VS-SEG-135`,
`VS-SEG-154`, `VS-SEG-188` (among the smallest patients, to fit the disk). Cite the dataset and paper if you use
results derived from them.

```bash
# dry run first: shows series, MB, licence and free disk; downloads nothing
python scripts/download_tcia.py Vestibular-Schwannoma-SEG data/raw/vs-seg VS-SEG-135 VS-SEG-154 VS-SEG-188
python scripts/download_tcia.py Pancreatic-CT-CBCT-SEG data/raw/pancreas-ct-cbct-seg \
    Pancreas-CT-CB_014 Pancreas-CT-CB_017
# add --yes to download. Each series is extracted atomically into its own directory and marked
# `.complete` only after the zip verified and the file count matched; an interrupted run simply
# resumes. A MANIFEST_<patient>.json is written once all of a patient's series succeeded.
```

### Citations and licence (required by TCIA's data usage policy)

Both datasets are licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/): you may share
and adapt them with attribution and an indication of changes. This repository does **not**
redistribute any of the data; figures and tables derived from it (metrics computed by this code,
plots without anatomy) are changed/derived material. If you use results derived from them, cite:

- Hong, J., Reyngold, M., Crane, C., Cuaron, J., Hajj, C., Mann, J., Zinovoy, M., Yorke, E.,
  LoCastro, E., Apte, A. P., & Mageras, G. (2021). *Breath-hold CT and cone-beam CT images with
  expert manual organ-at-risk segmentations from radiation treatments of locally advanced
  pancreatic cancer (Pancreatic-CT-CBCT-SEG)* (Version 2) [Data set]. The Cancer Imaging Archive.
  https://doi.org/10.7937/TCIA.ESHQ-4D90
- Shapey, J., Kujawa, A., Dorent, R., Wang, G., Bisdas, S., Dimitriadis, A., Grishchuk, D., Paddick, I.,
  Kitchen, N., Bradford, R., Saeed, S., Ourselin, S., & Vercauteren, T. (2021). *Segmentation of
  Vestibular Schwannoma from Magnetic Resonance Imaging: An Open Annotated Dataset and Baseline
  Algorithm (Vestibular-Schwannoma-SEG)* (Version 2) [Data set]. The Cancer Imaging Archive.
  https://doi.org/10.7937/TCIA.9YTJ-5Q73
- Clark, K., Vendt, B., Smith, K., et al. (2013). The Cancer Imaging Archive (TCIA): Maintaining and
  Operating a Public Information Repository. *Journal of Digital Imaging*, 26(6), 1045-1057.
  https://doi.org/10.1007/s10278-013-9622-7

These are also listed in [CITATION.cff](../CITATION.cff). Related papers: Han et al., *Med Phys* 2021
48(6):3084-3095 (pancreas methods paper); Shapey et al., *Sci Data* 2021 8:286 (VS-SEG).

### Using them

DICOM downloads do not say which structure set belongs to which image, or which dose to which plan.
`discover` follows the references stored in the files (RTDOSE → RTPLAN → RTSTRUCT → image series):

```bash
doseshift discover data/raw/vs-seg/VS-SEG-154 --roles TV=ptv,Cochlea=oar \
    --write-config experiments/real/vs_seg_154_audit.yaml
doseshift run experiments/real/vs_seg_154_audit.yaml
```

Real-data runs do not render the image underlay (`show_image: false` by default) so outputs do not
embed patient anatomy; set `show_image: true` in a config to see it locally. Ready-made configs for the five downloaded patients are in [`experiments/real/`](../experiments/real)
(`*_audit.yaml`: grid-assumption audit; `*_contours.yaml`: measured contour variation).

### Things to know about these data

- **Pairing matters.** Each VS-SEG patient has *two* RTDOSE and *two* RTSTRUCT (one set per MR
  series: contrast T1 and T2); picking by file name gives a silently wrong dose/contour pair. During
  development a wrong guess did exactly that, which is why `discover` and the "RTSTRUCT does not
  reference this image" warning exist.
- **No prescription in the pancreas data** (no RTPLAN) and no PTV structure. The configs use the
  median dose in the `ROI` structure (46 and 39 Gy) as a *reference dose*, labelled as such in
  every figure; Rx-dependent metrics (CI, GI, V%Rx) are omitted. VS-SEG's 13 Gy comes from the RTPLAN.
- **Pancreas CBCT-day contours** live in the planning frame of reference ("Aligned CT" series) and
  pair with the *planning* dose. That is what `contour_variation` uses. The original CBCT series
  (different frames of reference, 0.51 mm in-plane) are downloaded but not used yet.
- VS-SEG structures are tiny (tumour ~1-2 cc, cochlea 0.1 cc) on a 1 mm dose grid; contours include
  a `*Skull` ROI that is ignored.

## Not downloaded, and why

| Candidate | Status |
|---|---|
| HNSCC-3DCT-RT, MD Anderson HNSCC (TCIA) | **not in TCIA's public API listing**; may need a login or a signed data-use agreement. Accept those terms yourself on the TCIA site if you want them |
| HMM-RT / AIRTP (Hugging Face, `Jungle15/Radiotherapy_HaN_Lung_AIRTP`) | **gated**: manual approval after a form (real name, institutional email, affiliation), CC BY-NC-SA 4.0. Largest head-neck/lung CT+RTSTRUCT+RTPLAN+RTDOSE set; request access yourself |
| rest of the two collections above | disk space; the script fetches any further patient |

## No data needed

The analytic phantom drives every demo, test and the app. Export it as DICOM-RT to exercise the
reader (it writes CT, RTSTRUCT with series references, RTPLAN with the prescription, and RTDOSE):

```bash
doseshift make-phantom data/phantom                        # CT/, RTSTRUCT.dcm, RTPLAN.dcm, RTDOSE.dcm
doseshift make-phantom data/phantom_rev --reverse-dose-z   # RTDOSE frames in descending z
doseshift inspect --ct-dir data/phantom/CT --dose data/phantom/RTDOSE.dcm --rtstruct data/phantom/RTSTRUCT.dcm
```

## Reader requirements

One image series (CT or MR) with uniform slice spacing; absolute-Gy RTDOSE with uniform
`GridFrameOffsetVector` (either direction); planar RTSTRUCT contours. MR values are raw intensities
(no HU). See [docs/conventions.md](../docs/conventions.md).
