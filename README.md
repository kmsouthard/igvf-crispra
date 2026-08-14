# igvf-crispra

Element-level differential expression results for the **Hs27 fibroblast CRISPRa Perturb-seq**
screen, formatted for submission to the [IGVF Data Portal](https://data.igvf.org/).

This repository holds the code that converts the screen's mean-population matrix into the
element-level table IGVF expects, together with the promoter-window annotation it depends on and
the results of the run of record (2026-03-03).

Norman Lab, Memorial Sloan Kettering Cancer Center.

---

## What an "element" is here

An element is a **promoter window**, keyed on the tuple

```
(intended_target_name, intended_target_chr, intended_target_start, intended_target_end)
```

A gene with two annotated promoter windows therefore stays two distinct elements. In the run of
record, 307 elements map onto 291 distinct target genes.

Each element is represented by a **single guide** — the one with the highest `de_genes_fibro`
among that element's hit guides. The element-level table is that representative guide's effect
profile, not an aggregate across guides.

---

## Method

### Upstream (not in this repository)

Produced by **[norman-lab-msk/TFs_CRISPRa](https://github.com/norman-lab-msk/TFs_CRISPRa)** —
see `Code/Processing CRISPRa experiments and performing regressions/`, in particular *Step 3 —
Regressions and identifying masked_active guides for Hs27 experiment*. That repository, and the
paper's Methods section, are authoritative for everything in this subsection.

Libraries were aligned with **Cell Ranger v7.1.0** against **GRCh38-2020-A**, and guides called at a
5-UMI threshold.

**Normalization** (`offset_p_normalize`, restricted to stably captured genes, `mean > 0.2`) runs
per gem group:

1. Counts are converted to per-cell transcript probabilities by dividing by each cell's UMI total.
2. Droplet occupancy is corrected by stratifying cells on `number_of_cells` per droplet and
   computing the average probability per occupancy class using a **2 %-trimmed mean**, then dividing
   each cell by its class scale factor.
3. Values are rescaled by `pairwise_singlet_equivalent_UMI_count` to adjust transcript capture.
4. The **mean of control cells in the same gem group** is subtracted.
5. The result is divided by the **standard deviation of control cells in the same gem group**.

Centering and scaling are therefore both done against controls *within each gem group*, not against
a single pooled control distribution. This is what puts effect scores in units of control standard
deviation before any regression is fitted.

**Regression** is ordinary least squares (`scipy.linalg.lstsq`) of that standardized matrix on the
guide design matrix. The reported effect score is the OLS **coefficient**. P-values come from a
two-sided *t*-test on `coef / SE` with `df = n − rank`, and Benjamini–Hochberg correction across
genes produces the adjusted values.

**Active-guide calls** are made on **target-masked** coefficients: each guide's own target gene
column is set to `NaN` so that on-target activation cannot drive the call. Guides are then clustered
on `1 − pairwise correlation` of those masked profiles with average-linkage hierarchical clustering,
cut by the inconsistency criterion. A guide is `masked_active` when its cluster contains **two or
more sgRNAs against the same gene** — reproducibility across independent guides, not effect size.
`expanded_masked_active` is a disjoint rescue set from revisiting the clustering at threshold 1.2,
keeping good, non-sequence-driven guides whose target gene no cluster had yet recovered.

The result is `fibroblast_CRISPRa_mean_pop.h5ad` (10,916 guides × 4,914 genes), whose `X` holds
effect scores and whose `p`, `adj_p` and `masked` layers hold p-values, adjusted p-values and the
target-masked coefficients.

> **This pipeline reads `X`, not `masked`.** Effect scores in the output are the unmasked
> coefficients, so a guide's own target gene appears with its real on-target effect. The masking
> exists only to make the *activity classification* independent of on-target signal. Comparing these
> effect scores to the clustering in the paper without knowing that will mislead.

> Note the Cell Ranger version: this screen used **7.1.0 / GRCh38-2020-A**, which differs from the
> lab's dual-guide fibroblast work (9.0.1 / GRCh38-2024-A). Do not carry one version across to the
> other when registering analysis steps.

### This repository

`generate_element_level_results.py` performs four steps.

1. **Join guide activity metadata.** Table S4 (`Guide Activity`) is joined onto `adata.obs` by
   `guide_identity`, supplying `seed_driven_fibro`, `bad_seed`, `active_fibro`,
   `expanded_active_fibro` and `de_genes_fibro`.

   These are the published names for the upstream flags, and they do **not** all correspond
   one-to-one with the `obs` columns already in the h5ad. Verified against this h5ad:

   | Table S4 column | h5ad `obs` column | Relationship |
   |---|---|---|
   | `active_fibro` (659 true) | `masked_active` (661 true) | identical on all 10,707 joined guides |
   | `expanded_active_fibro` (106 true) | not in h5ad | disjoint from `active_fibro`, zero overlap |
   | `seed_driven_fibro` (1,038 true) | `sequence_driven` (1,083 true) | **differ on 454 guides — not interchangeable** |

   Table S4 covers 10,707 of the 10,916 guides in the h5ad. The remaining 209 are outside the
   classified set by design — 17 `off-target` controls, plus 192 guides (six each) against 32 TF
   genes that were not carried into the activity table. They join as `NaN` and are excluded by the
   hit filter, which is the intended behaviour.

2. **Select hit guides.**

   ```
   seed_driven_fibro == False
   AND bad_seed == False
   AND (active_fibro OR expanded_active_fibro)
   ```

3. **Pick one representative guide per element.** Hit guides are joined to their promoter windows,
   grouped by the element tuple, and the guide with the highest `de_genes_fibro` is kept.

4. **Emit long-form results.** Effect scores and both p-value layers are melted to one row per
   element × gene, gene symbols are mapped to unversioned Ensembl IDs, and rows are filtered to
   `p_val_adj < 0.05`.

### Cohort sizes, run of record

| Stage | Count |
|---|---|
| Guides in the mean-pop h5ad | 10,916 |
| Classified in Table S4 (rest are controls / unclassified targets) | 10,707 |
| `active_fibro` / `expanded_active_fibro` (disjoint) | 659 / 106 |
| Hit guides after seed / bad-seed filter | 626 |
| Hit guides matched to a promoter window | 626 (none lost) |
| Distinct promoter elements | 313 |
| Elements present in output | 307 |
| Distinct target genes | 291 |
| Element × gene rows written | 31,047 |

---

## Inputs

Two inputs are shipped here; two are too large or not ours to redistribute. All four are pinned by
SHA-256 so a future run can prove it used the same data.

| File | Shipped | SHA-256 |
|---|---|---|
| `merged_guides_promoters_UPSTREAM_PRIORITIZED.csv` | yes | `1ac915099d18032148841b862331ce7bb33d24bf811f9bee01424c4b534c6f0f` |
| `Supplementary_Table1.xlsx` | no | `dc125f4727d650163fbfb08175d2a48ea77406787e77d3fa21a134639266a6c0` |
| `fibroblast_CRISPRa_mean_pop.h5ad` | no | `6e10a2605a3dbc448ac4756399adc947a493a75f8961d3793af8f851aa0411cf` |
| `hit_guides_by_genes.tsv` | no | not an input to this script |

**`Supplementary_Table1.xlsx`** (17 MB) is manuscript supplementary material and is distributed
with the paper rather than from here. Only the sheet `Table S4 Guide Activity` is read.

**`fibroblast_CRISPRa_mean_pop.h5ad`** (1.7 GB) lives at
`/data1/normantm/datasets/Hs27_CRISPRa/fibroblast_CRISPRa_mean_pop.h5ad` on the MSK cluster.

**`hit_guides_by_genes.tsv`** (275 MB) is a wide guide × gene matrix described in the file format
spec. It is a sibling deliverable, not an input to this script, and exceeds GitHub's 100 MB limit.

---

## Outputs

`results/element_level_results.tsv` — one row per element × gene, restricted to `p_val_adj < 0.05`.

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

The two `target_gene` / `intended_target_name` columns are easy to confuse: the first is the
readout, the second is the perturbation.

---

## Running it

```bash
python generate_element_level_results.py \
  --h5ad         /data1/normantm/datasets/Hs27_CRISPRa/fibroblast_CRISPRa_mean_pop.h5ad \
  --guide_meta   Supplementary_Table1.xlsx \
  --guide_coords merged_guides_promoters_UPSTREAM_PRIORITIZED.csv \
  --padj_thresh  0.05 \
  --out          results/element_level_results.tsv
```

The guide summary is written alongside `--out` with `_guide_summary` appended.

### Environment

Verified under `milo-env`: Python 3.13.7, scanpy 1.11.4, anndata 0.12.2, pandas 2.3.2, numpy 2.3.3.

The two notebooks record kernel `scanpy_perturbseq`, which currently fails to import pandas on
isce001; treat the versions above as the reproducible environment.

---

## Known behaviour worth knowing before you cite these files

**Elements with no significant gene are absent, not zero.** The `p_val_adj < 0.05` filter is applied
before the summary is built, so an element whose representative guide has no gene under threshold
drops out of both output files rather than appearing with `n_de_genes = 0`. In the run of record
that is 6 elements: 313 go in, 307 come out. A consumer of these files cannot distinguish a
tested-but-null element from one that was never tested.

This is deliberate for a hit-list deliverable, but if IGVF wants the full tested set, the filter
needs to move after the summary step.

**Two printed diagnostics cannot fire.** For the same reason, the console lines
`Elements with ≥1 DE gene` and `Elements with 0 DE genes` always report the full count and zero
respectively. They do not measure what their wording suggests. The written files are unaffected.

**Guide coordinates are computed but not emitted.** `Start` and `End` from the coordinates file are
renamed to `guide_start` / `guide_end` and then dropped before writing, since IGVF's element-level
layout keys on the promoter window rather than the protospacer.

---

## IGVF submission mapping

The output is a `tabular_file`. The content type depends on whether it is being submitted as the
element-level quantification or as a per-element differential expression product:

```
content_type: "differential element quantifications"
          or  "local differential expression per element"
```

Registering the `analysis_step_version` requires `lab`, `award`, `analysis_step` and at least one
`software_versions` entry. Point the `Software.source_url` at this repository and the
`SoftwareVersion` at the release tag.

---

## Contents

```
generate_element_level_results.py             the element-level generator
merged_guides_promoters_UPSTREAM_PRIORITIZED.csv   promoter windows per guide
conversion_of_h5ad_to_tsv.ipynb               h5ad → wide TSV conversion
tts_guide_merge.ipynb                         TSS / promoter annotation build
results/                                      run of record, 2026-03-03
```
