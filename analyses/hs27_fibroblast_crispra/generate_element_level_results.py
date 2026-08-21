import pandas as pd
import scanpy as sc
import argparse

def parse_args():
    p = argparse.ArgumentParser()
    # Download from Zenodo: https://doi.org/10.5281/zenodo.15200179
    p.add_argument('--h5ad', default='fibroblast_CRISPRa_mean_pop.h5ad')
    p.add_argument('--guide_meta', default='Supplementary_Table1.xlsx')
    p.add_argument('--guide_coords', default='merged_guides_promoters_UPSTREAM_PRIORITIZED.csv')
    p.add_argument('--padj_thresh', type=float, default=0.05)
    p.add_argument('--out', default='element_level_results.tsv')
    return p.parse_args()

def main():
    args = parse_args()

    adata = sc.read_h5ad(args.h5ad)

    # Merge guide activity metadata from Excel (not present in base h5ad)
    guide_meta = pd.read_excel(args.guide_meta, sheet_name='Table S4 Guide Activity')
    meta_cols = ['guide_identity', 'seed_driven_fibro', 'bad_seed',
                 'active_fibro', 'expanded_active_fibro', 'de_genes_fibro']
    guide_meta = guide_meta[[c for c in meta_cols if c in guide_meta.columns]].set_index('guide_identity')
    adata.obs = adata.obs.join(guide_meta, how='left')

    # Hit guide filter
    hit_mask = (
        (adata.obs.seed_driven_fibro == False) &
        (adata.obs.bad_seed == False) &
        (adata.obs.active_fibro | adata.obs.expanded_active_fibro)
    )
    hits_obs = adata.obs[hit_mask].copy()
    hits_idx = hits_obs.index

    # Merge promoter coordinates onto guide metadata
    coords = pd.read_csv(
        args.guide_coords,
        usecols=['guide_identity', 'ensembl', 'Chromosome', 'Start', 'End', 'Strand',
                 'protospacer', 'promoter_start', 'promoter_end']
    ).drop_duplicates('guide_identity').rename(columns={
        'guide_identity':  'guide_id',
        'ensembl':         'intended_target_name',
        'Chromosome':      'intended_target_chr',
        'promoter_start':  'intended_target_start',
        'promoter_end':    'intended_target_end',
        'Chromosome':      'intended_target_chr',
        'Start':           'guide_start',
        'End':             'guide_end',
        'Strand':          'strand',
        'protospacer':     'spacer',
    })

    hits_obs = hits_obs.merge(coords, left_index=True, right_on='guide_id', how='left')

    # Select one representative guide per promoter element:
    # the guide with the highest de_genes_fibro (most DE genes driven).
    # A promoter element is uniquely defined by (intended_target_name, intended_target_chr,
    # intended_target_start, intended_target_end) — handles genes with multiple promoter windows.
    element_cols = ['intended_target_name', 'intended_target_chr',
                    'intended_target_start', 'intended_target_end']
    rep_guides = (
        hits_obs.sort_values('de_genes_fibro', ascending=False)
                .groupby(element_cols, as_index=False)
                .first()  # top guide per element after sorting by de_genes_fibro
    )[['guide_id'] + element_cols]

    # Build long-form effect/p-value table for representative guides only
    rep_idx = rep_guides['guide_id'].values
    sym_to_ensg = adata.var['gene_id'].str.split('.').str[0].to_dict()

    def melt_layer(layer_key, value_name):
        df = adata[rep_idx].to_df(layer=layer_key) if layer_key else adata[rep_idx].to_df()
        return df.reset_index().melt(id_vars='guide_identity', var_name='gene_symbol', value_name=value_name)

    long = (
        melt_layer(None,    'effect_score')
        .merge(melt_layer('p',     'p_val'),     on=['guide_identity', 'gene_symbol'])
        .merge(melt_layer('adj_p', 'p_val_adj'), on=['guide_identity', 'gene_symbol'])
    )
    long['target_gene'] = long['gene_symbol'].map(sym_to_ensg)
    long = long.drop(columns='gene_symbol')
    long = long[long.p_val_adj < args.padj_thresh]
    long.rename(columns={'guide_identity': 'guide_id'}, inplace=True)

    # Merge element metadata back onto long-form results
    long = long.merge(rep_guides, on='guide_id', how='left')

    for col in ['intended_target_start', 'intended_target_end']:
        long[col] = pd.to_numeric(long[col], errors='coerce').astype('Int64')

    out_cols = ['effect_score', 'p_val', 'p_val_adj', 'guide_id', 'target_gene',
                'intended_target_name', 'intended_target_chr',
                'intended_target_start', 'intended_target_end']
    long[out_cols].to_csv(args.out, sep='\t', index=False)
    print(f"Wrote {len(long)} element×gene pairs ({long['guide_id'].nunique()} elements) to {args.out}")

    # Summary: DE gene counts per representative guide
    summary = (
        long.groupby(['guide_id', 'intended_target_name'])
            .agg(n_de_genes=('target_gene', 'count'),
                 n_upregulated=('effect_score', lambda x: (x > 0).sum()),
                 n_downregulated=('effect_score', lambda x: (x < 0).sum()))
            .reset_index()
            .merge(rep_guides[['guide_id', 'intended_target_chr',
                                'intended_target_start', 'intended_target_end']], on='guide_id')
            .sort_values('n_de_genes', ascending=False)
    )
    summary_path = args.out.replace('.tsv', '_guide_summary.tsv')
    summary.to_csv(summary_path, sep='\t', index=False)

    print(f"\nDE gene summary (p_val_adj < {args.padj_thresh}):")
    print(f"  Elements with ≥1 DE gene:  {(summary.n_de_genes >= 1).sum()} / {len(summary)}")
    print(f"  Median DE genes/element:   {summary.n_de_genes.median():.0f}")
    print(f"  Mean DE genes/element:     {summary.n_de_genes.mean():.1f}")
    print(f"  Max DE genes/element:      {summary.n_de_genes.max()} ({summary.iloc[0].guide_id})")
    print(f"  Elements with 0 DE genes:  {(summary.n_de_genes == 0).sum()}")
    print(f"\nTop 10 elements by DE gene count:")
    print(summary[['guide_id', 'intended_target_name', 'n_de_genes',
                   'n_upregulated', 'n_downregulated']].head(10).to_string(index=False))
    print(f"\nSummary written to {summary_path}")

if __name__ == '__main__':
    main()