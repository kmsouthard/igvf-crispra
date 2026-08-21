# RPE-1 Multiome Perturb-seq (CRISPRi) — analysis code

IGVF-facing analysis code for the **hTERT RPE-1 Multiome Perturb-seq** screen.

**Paper.** Metzner E\*, Southard KM\*, Norman TM. *Multiome Perturb-seq unlocks scalable discovery
of integrated perturbation effects on the transcriptome and epigenome.* Cell Systems **16**, 101161
(2025). [10.1016/j.cels.2024.12.002](https://doi.org/10.1016/j.cels.2024.12.002) (\*equal
contribution)

> **Status: scaffold.** No code added yet.

## The screen

CRISPRi (dCas9-ZIM3) in hTERT RPE-1, reading out gene expression and chromatin accessibility from
the same nuclei via 10x Multiome. **13 chromatin remodelers** — SMARCE1, SMARCB1, ARID1A, SMARCA4,
DPF2, SMARCC1, SMARCC2, EP400, ACTL6A, DMAP1, SUZ12, EZH2, YY1 — plus 3 non-targeting controls.

**Open question, blocks everything downstream:** how an IGVF "element" is defined on the ATAC side —
peak, or perturbed promoter. Unlike the [Hs27 CRISPRa analysis](../hs27_fibroblast_crispra/), this
screen has two modalities, and the answer sets the output schema. Document it here once settled.

## Data availability

| Resource | Identifier |
|---|---|
| Raw sequencing | SRA [PRJNA1128171](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1128171) |
| Paper's full analysis code | [10.5281/zenodo.14217682](https://doi.org/10.5281/zenodo.14217682) |

IGVF portal objects — GEX measurement set `IGVFDS4073QPQB`, ATAC measurement set `IGVFDS9292PROA`,
guide sequencing auxiliary set `IGVFDS4833CNYU`, guide RNA sequences file `IGVFFI1270VQEU`, curated
set `IGVFDS2221DOTB`.

## What belongs here

Only what the IGVF submission depends on — the code that produces the portal tables, plus its
annotation inputs — so the `AnalysisStepVersion` has something citable to point at. The paper's full
analysis is at the Zenodo DOI above and doesn't need duplicating.

Layout mirrors the directory next door: `notebooks/`, `results/`, and this README carrying the
method, SHA-256-pinned inputs and output schema. Repo-wide conventions are in the
[root README](../../README.md).
