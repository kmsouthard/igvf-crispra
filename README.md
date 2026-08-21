# igvf-perturbseq

Analysis code behind the Norman Lab's **Perturb-seq** submissions to the
[IGVF Data Portal](https://data.igvf.org/)

| Directory | Screen | Paper |
|---|---|---|
| [`analyses/hs27_fibroblast_crispra/`](analyses/hs27_fibroblast_crispra/) | Hs27 fibroblast, 1,836 TFs, CRISPRa (RNA) | Southard and Ardy et al., *Nat Genet* **57**, 2323–2334 (2025), [10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1) |
| [`analyses/rpe1_multiome_crispri/`](analyses/rpe1_multiome_crispri/) | hTERT RPE-1, 13 chromatin remodelers, CRISPRi (RNA + ATAC) | Metzner and Southard et al., *Cell Syst* **16**, 101161 (2025), [10.1016/j.cels.2024.12.002](https://doi.org/10.1016/j.cels.2024.12.002) |

## Environment

`environment.yml` / `requirements.txt` are the shared baseline.

```bash
conda env create -f environment.yml      # or: pip install -r requirements.txt
```
