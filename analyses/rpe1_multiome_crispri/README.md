# RPE-1 Multiome Perturb-seq (CRISPRi) — analysis code

**Paper.** Metzner E*, Southard KM*, Norman TM (*equal contribution).
*Multiome Perturb-seq unlocks scalable discovery of integrated perturbation effects on the
transcriptome and epigenome.* Cell Systems **16**, 101161 (2025).
[10.1016/j.cels.2024.12.002](https://doi.org/10.1016/j.cels.2024.12.002)

> **Scaffold — no code added yet.** Only what the IGVF submission depends on belongs here; the
> paper's full analysis is at the Zenodo DOI below.

## Screen

CRISPRi (dCas9-ZIM3) in hTERT RPE-1, gene expression and chromatin accessibility from the same
nuclei via 10x Multiome. 13 chromatin remodelers — SMARCE1, SMARCB1, ARID1A, SMARCA4, DPF2,
SMARCC1, SMARCC2, EP400, ACTL6A, DMAP1, SUZ12, EZH2, YY1 — plus 3 non-targeting controls.

## "element" definition

**Unsettled — peak, or perturbed promoter.** Two modalities, unlike the
[Hs27 CRISPRa analysis](../hs27_fibroblast_crispra/) where an element is a promoter window. Sets
the output schema, so it blocks results.

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
