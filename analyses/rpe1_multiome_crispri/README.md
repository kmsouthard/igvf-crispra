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

Run of record, 2026-08-21, in `results/`.

| File | Grain | Rows |
|---|---|---|
| `multiome_paper_guide_effect_matrix.csv` | element × gene | 14,868 |
| `multiome_paper_differential_peaks_by_guide.csv` | element × peak | 10,935 |

Columns follow the Hs27 schema — `effect_score`, `p_val`, `p_val_adj`, `guide_id`,
`intended_target_name`, `intended_target_chr/_start/_end` — plus `target_gene` (GEX) or
`chr/start/end` of the peak (ATAC).

## Environment

`environment.yml` in this directory, on top of the repository baseline. **snapatac2 must be
2.6.0** — other versions do not reproduce the published Table S3. `perturbseq` is not packaged;
see [thomasmaxwellnorman/perturbseq_demo](https://github.com/thomasmaxwellnorman/perturbseq_demo).

## Inputs

Pinned by SHA-256.

| File | Shipped | SHA-256 |
|---|---|---|
| `notebooks/multiome_paper_igvf_guides.csv` — 16 guides, the seed table | yes | `beaa98924e3f78a5b2ac12bce14a5e64e746c9e7f0e1d3652b7c68cec81deb41` |
| `atac_singlets_macs3_peaks.h5ad` — 1.6 GB | no | `90dbc9f9613d09a146cbd60d4376a84debc1cf05b0478456fe750b9c265fd45c` |
| `gex_norm_regressed.hdf5` — 370 MB | no | `2ba7f3ac585f8e2b04be423de8a8c2fd3a38a441c3f769a80061167f2865bd71` |
| `epdNewHuman006_extended_promoter_regions.bed` — 1.7 MB | no | `b6777b8e0b78b8cc06e1f9fe710eba8e5feae661beb30904a14040832443d300` |
| `genes.gtf` — 1.4 GB | no | `2dc6e7406e883a146c7cc933a2b08c8d0546e7b57e0487a93cbbc1c455868528` |
| `genome.fa` — 3.0 GB | no | `fb7421217e7058120cd60a5277445198e8deef2dff1edc46cd1e98b31fe64cbb` |
| `perturbseq` package — 7 `.py` files, content digest | no | `dfeb9ee6dc4c04c5bbfff4579a054e25f6abd082101266305ca968eb1a8552d2` |

**Data availability.** Raw sequencing: SRA
[PRJNA1128171](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1128171). Paper's full analysis
code: [10.5281/zenodo.14217682](https://doi.org/10.5281/zenodo.14217682).

## IGVF portal objects

| Object | Accession |
|---|---|
| GEX measurement set | `IGVFDS4073QPQB` |
| ATAC measurement set | `IGVFDS9292PROA` |
| Guide sequencing auxiliary set | `IGVFDS4833CNYU` |
| Guide RNA sequences file | `IGVFFI1270VQEU` |
| Curated set | `IGVFDS2221DOTB` |
