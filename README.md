# igvf-crispra

Element-level differential expression results for the **Hs27 fibroblast CRISPRa Perturb-seq**
screen, as submitted to the [IGVF Data Portal](https://data.igvf.org/).

Holds the code that converts the screen's mean-population matrix into IGVF's element-level table,
the promoter-window annotation it depends on, and the results of record (2026-03-03).

Norman Lab, Memorial Sloan Kettering Cancer Center.

---

## What an "element" is here

A **promoter window**, keyed on
`(intended_target_name, intended_target_chr, intended_target_start, intended_target_end)`. A gene
with two annotated windows stays two elements — 307 elements map onto 291 genes.

Each element is represented by a **single guide**, the one with the highest `de_genes_fibro` among
its hit guides. The table is that guide's effect profile, not an aggregate.

---

## Method

### Upstream (not in this repository)

Produced by **[norman-lab-msk/TFs_CRISPRa](https://github.com/norman-lab-msk/TFs_CRISPRa)**, in
`Code/Processing CRISPRa experiments and performing regressions/` — chiefly *Step 3 — Regressions
and identifying masked_active guides for Hs27 experiment*. That repository and the paper's Methods
are authoritative here; the summary below is orientation only.

Cell Ranger **v7.1.0** against **GRCh38-2020-A**, guides called at a 5-UMI threshold.

**Normalization** (`offset_p_normalize`, on stably captured genes at `mean > 0.2`) runs per gem
group: counts become per-cell transcript probabilities; droplet occupancy is corrected by
stratifying on cells-per-droplet and dividing by a **2 %-trimmed mean** per occupancy class; values
are rescaled by `pairwise_singlet_equivalent_UMI_count`; then the **mean of control cells in the
same gem group** is subtracted and the result divided by **those controls' standard deviation**.

Centering and scaling are both against controls *within each gem group*, not a pooled control
distribution — that is what puts effect scores in control-SD units before any regression is fitted.

**Regression** is ordinary least squares (`scipy.linalg.lstsq`) of that matrix on the guide design
matrix. The effect score is the OLS coefficient; p-values are a two-sided *t*-test on `coef / SE`
with `df = n − rank`, Benjamini–Hochberg corrected across genes.

**Active-guide calls** use **target-masked** coefficients — each guide's own target column set to
`NaN`, so on-target activation cannot drive the call. Guides are clustered on `1 − correlation` of
those masked profiles (average linkage, inconsistency criterion); a guide is `masked_active` when
its cluster holds **two or more sgRNAs against the same gene**. That is a reproducibility criterion,
not an effect-size one. `expanded_masked_active` is a disjoint rescue set from re-clustering at
threshold 1.2.

The result is `fibroblast_CRISPRa_mean_pop.h5ad` (10,916 guides × 4,914 genes): `X` holds effect
scores, layers `p`, `adj_p` and `masked` hold p-values, adjusted p-values and the masked
coefficients.

> **This pipeline reads `X`, not `masked`** — so a guide's own target gene carries its real
> on-target effect. Masking exists only to keep the *activity classification* independent of
> on-target signal.

### This repository

`generate_element_level_results.py`:

1. **Join guide activity metadata** — Table S4 (`Guide Activity`) onto `adata.obs` by
   `guide_identity`, supplying `seed_driven_fibro`, `bad_seed`, `active_fibro`,
   `expanded_active_fibro`, `de_genes_fibro`.
2. **Select hit guides** — `seed_driven_fibro == False AND bad_seed == False AND (active_fibro OR
   expanded_active_fibro)`.
3. **Pick one representative guide per element** — highest `de_genes_fibro` within each element.
4. **Emit long-form results** — effect scores and both p-value layers melted to one row per
   element × gene, symbols mapped to unversioned Ensembl IDs, filtered to `p_val_adj < 0.05`.

### Cohort sizes, run of record

| Stage | Count |
|---|---|
| Guides in the mean-pop h5ad | 10,916 |
| Classified in Table S4 | 10,707 |
| `active_fibro` / `expanded_active_fibro` (disjoint) | 659 / 106 |
| Hit guides after seed / bad-seed filter | 626 |
| Hit guides matched to a promoter window | 626 (none lost) |
| Distinct promoter elements | 313 |
| Elements present in output | 307 |
| Distinct target genes | 291 |
| Element × gene rows written | 31,047 |

---

## Inputs

Pinned by SHA-256 so a future run can prove it used the same data.

| File | Shipped | SHA-256 |
|---|---|---|
| `merged_guides_promoters_UPSTREAM_PRIORITIZED.csv` | yes | `1ac915099d18032148841b862331ce7bb33d24bf811f9bee01424c4b534c6f0f` |
| `Supplementary_Table1.xlsx` — 17 MB, manuscript material, distributed with the paper. Only sheet `Table S4 Guide Activity` is read. | no | `dc125f4727d650163fbfb08175d2a48ea77406787e77d3fa21a134639266a6c0` |
| `fibroblast_CRISPRa_mean_pop.h5ad` — 1.7 GB, from Zenodo (see below). | no | `6e10a2605a3dbc448ac4756399adc947a493a75f8961d3793af8f851aa0411cf` |

`hit_guides_by_genes.tsv` (275 MB) is a sibling deliverable described in the file format spec, not
an input to this script, and exceeds GitHub's 100 MB limit.

### Data availability

The mean-population matrix is published in the **Hs27-CRISPRa-TFs** record,
[10.5281/zenodo.15200179](https://doi.org/10.5281/zenodo.15200179) — the copy there is byte-identical
to the one used for the run of record (md5 `ba44c7813903bb5df900348d6b0d589a`).

Other processed datasets for this screen are in the lab's Zenodo community,
<https://zenodo.org/communities/normanlabmsk/>; raw FASTQs are at SRA
[PRJNA1108254](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1108254).

---

## Outputs

`results/element_level_results.tsv` — one row per element × gene at `p_val_adj < 0.05`.

| Column | Meaning |
|---|---|
| `effect_score` | Deviation from control in SD units. Negative is down. |
| `p_val` | Uncorrected p-value from the regression. |
| `p_val_adj` | Benjamini–Hochberg adjusted p-value. |
| `guide_id` | Representative guide for the element, `SYMBOL_PROTOSPACER`. |
| `target_gene` | Ensembl ID of the **measured** gene whose expression changed. |
| `intended_target_name` | Ensembl ID of the gene whose promoter was **perturbed**. |
| `intended_target_chr` | Chromosome of the promoter window, UCSC style (`chr4`). |
| `intended_target_start` | Promoter window start. |
| `intended_target_end` | Promoter window end. |

`results/element_level_results_guide_summary.tsv` — one row per element, with `n_de_genes`,
`n_upregulated`, `n_downregulated` and the element coordinates.

`target_gene` is the readout and `intended_target_name` is the perturbation — easy to transpose.

---

## Running it

Fetch the mean-population matrix from Zenodo, then run from the repository root — every other input
and every default path is relative to it.

```bash
conda env create -f environment.yml      # or: pip install -r requirements.txt

python generate_element_level_results.py \
  --h5ad         fibroblast_CRISPRa_mean_pop.h5ad \
  --guide_meta   Supplementary_Table1.xlsx \
  --guide_coords merged_guides_promoters_UPSTREAM_PRIORITIZED.csv \
  --padj_thresh  0.05 \
  --out          results/element_level_results.tsv
```

Those are also the defaults, so a bare `python generate_element_level_results.py` does the same
thing. The guide summary is written alongside `--out` with `_guide_summary` appended.

Tested with Python 3.13.7, scanpy 1.11.4, anndata 0.12.2, pandas 2.3.2, numpy 2.3.3, scipy 1.16.2,
h5py 3.14.0, openpyxl 3.1.5.

---

## Contents

```
generate_element_level_results.py                    the element-level generator
merged_guides_promoters_UPSTREAM_PRIORITIZED.csv     promoter windows per guide
conversion_of_h5ad_to_tsv.ipynb                      h5ad → wide TSV conversion
tts_guide_merge.ipynb                                TSS / promoter annotation build
results/                                             run of record, 2026-03-03
environment.yml / requirements.txt                   conda and pip dependencies
```
