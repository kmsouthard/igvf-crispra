# igvf-perturbseq

Analysis code behind the Norman Lab's **Perturb-seq** submissions to the
[IGVF Data Portal](https://data.igvf.org/) — what turns each screen's processed matrices into the
tables IGVF asks for, plus the results of record.

Norman Lab, Memorial Sloan Kettering Cancer Center.

## Analyses

Each directory is self-contained — method, inputs, outputs and how to run it are in its own README.
The two screens define an "element" differently and are not interchangeable.

| Directory | Screen | Paper |
|---|---|---|
| [`analyses/hs27_fibroblast_crispra/`](analyses/hs27_fibroblast_crispra/) | Hs27 fibroblast, 291 TFs, CRISPRa (RNA) | Southard et al., *Nat Genet* **57**, 2323–2334 (2025), [10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1) |
| [`analyses/rpe1_multiome_crispri/`](analyses/rpe1_multiome_crispri/) | hTERT RPE-1, 13 chromatin remodelers, CRISPRi (RNA + ATAC) | Metzner et al., *Cell Syst* **16**, 101161 (2025), [10.1016/j.cels.2024.12.002](https://doi.org/10.1016/j.cels.2024.12.002) |

## Environment

`environment.yml` / `requirements.txt` are the shared baseline; an analysis needing more ships its
own and says so in its README.

```bash
conda env create -f environment.yml      # or: pip install -r requirements.txt
```

## Conventions

The IGVF portal hangs a `Software` record and an `AnalysisStepVersion` off a URL and a tag here, so
history has to stay stable.

- **Branch and PR.** Never force-push `main` — a rewritten tag breaks a portal record.
- **Tag every submitted run**, and name the tag in the analysis README.
- **One directory per screen**, owned by whoever ran it. Shared code moves to a top-level module
  only once a second analysis actually uses it.
- **Pin inputs by SHA-256** in the analysis README, including inputs too large to commit.
- **`results/` is a dated run of record** — regenerating means a new commit, not a silent overwrite.
- **Never commit matrices.** `.h5ad`, `.h5`, `.hdf`, `.loom`, `.mtx` are ignored repo-wide; anything
  over ~50 MB needs a deposit and a checksum instead.
