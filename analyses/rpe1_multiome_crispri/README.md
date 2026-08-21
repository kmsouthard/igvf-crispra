# RPE-1 Multiome Perturb-seq (CRISPRi) — analysis code

**Paper.** Metzner E*, Southard KM*, Norman TM (*equal contribution).
*Multiome Perturb-seq unlocks scalable discovery of integrated perturbation effects on the
transcriptome and epigenome.* Cell Systems **16**, 101161 (2025).
[10.1016/j.cels.2024.12.002](https://doi.org/10.1016/j.cels.2024.12.002)

> Only what the IGVF submission depends on belongs here; the paper's full analysis is at the
> Zenodo DOI below.

## Screen

CRISPRi (dCas9-ZIM3) in hTERT RPE-1, gene expression and chromatin accessibility from the same
nuclei via 10x Multiome. 13 chromatin remodelers — SMARCE1, SMARCB1, ARID1A, SMARCA4, DPF2,
SMARCC1, SMARCC2, EP400, ACTL6A, DMAP1, SUZ12, EZH2, YY1 — plus 3 non-targeting controls.

## "element" definition

The **perturbed promoter**, keyed on `intended_target_chr/_start/_end` — nearest EPD extended
promoter window to the guide, same convention as the
[Hs27 CRISPRa analysis](../hs27_fibroblast_crispra/). The two modalities differ in the *readout*,
not the element: GEX reports `target_gene`, ATAC reports the differential peak's `chr/start/end`.

## Processing

`notebooks/multiome_data_organization.ipynb`:

1. Maps each guide protospacer to genomic coordinates, then to its nearest EPD promoter window.
2. **GEX** — `ks_de` against NTC cells on the normalized/regressed population; effect score is the
   mean-population z-score, filtered to `q < 0.1`, symbols mapped to Ensembl IDs via the GTF.
3. **ATAC** — MACS3 peaks merged, paired-insertion peak matrix, Mann–Whitney U per guide against
   NTC with Benjamini–Hochberg control; effect score is log2 fold-change, filtered to `q < 0.1`.

Note both modalities use `q < 0.1`, where the Hs27 analysis uses `p_val_adj < 0.05`.

## Outputs

| File | Grain |
|---|---|
| `multiome_paper_guide_effect_matrix.csv` | element × gene |
| `multiome_paper_differential_peaks_by_guide.csv` | element × peak |

Columns follow the Hs27 schema — `effect_score`, `p_val`, `p_val_adj`, `guide_id`,
`intended_target_name`, `intended_target_chr/_start/_end` — plus `target_gene` (GEX) or
`chr/start/end` of the peak (ATAC).

## Environment

`environment.yml` in this directory, on top of the repository baseline. **snapatac2 must be
2.6.0** — other versions do not reproduce the published Table S3. `perturbseq` is not packaged;
see [thomasmaxwellnorman/perturbseq_demo](https://github.com/thomasmaxwellnorman/perturbseq_demo).

## Inputs

| Resource | Identifier |
|---|---|
| Raw sequencing | SRA [PRJNA1128171](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1128171) |
| Paper's analysis code | [10.5281/zenodo.14217682](https://doi.org/10.5281/zenodo.14217682) |

## IGVF portal objects

| Object | Accession |
|---|---|
| GEX measurement set | `IGVFDS4073QPQB` |
| ATAC measurement set | `IGVFDS9292PROA` |
| Guide sequencing auxiliary set | `IGVFDS4833CNYU` |
| Guide RNA sequences file | `IGVFFI1270VQEU` |
| Curated set | `IGVFDS2221DOTB` |
