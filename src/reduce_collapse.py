#!/usr/bin/env python3
# python/collapse.py

import pandas as pd
import numpy as np
import argparse

def load_data(guides_file):
    df = pd.read_csv(guides_file)
    return df

def sum_expression(x):
    # 1. Handle the "All Missing" case
    # If all values are -1, the entire guide targets unmeasured transcripts
    if len(x) > 0 and (x < 0).all():
        return -1.0000
    
    # 2. Handle the "Mixed" or "Valid" case
    # Filter for values >= 0 (0.0 to 1.0) to ignore the -1 flags
    valid_values = x[x >= 0]
    
    # If there are valid values, sum only those. 
    # If no valid values exist (but the list wasn't all < 0), return -1
    if len(valid_values) > 0:
        return round(valid_values.sum(), 4)
    
    return -1.0000

def collapse_guides(expanded_df):
    # Ensure median_expression_norm is present
    if 'median_expression_norm' not in expanded_df.columns:
        expanded_df['median_expression_norm'] = -1.0000

    # Sort by Guide Sequence and transcript_id to ensure first alphabetically
    sorted_df = expanded_df.sort_values(['Guide Sequence', 'Symbol'])
    sorted_df['Symbol_Contig_Idx'] = sorted_df['Symbol'].astype(str) + ':' + sorted_df['Contig_Idx'].astype(str)
    
    collapse_df = sorted_df.groupby('Guide Sequence', as_index=False).agg(
        **{
            'Target Sequence': ('Target Sequence', 'first'),
            'Guide Score': ('Guide Score', 'max'),
            'Title': ('Title', 'first'),
            'Gene': ('Gene', 'first'),
            'Symbol': ('Symbol', lambda x: '|'.join(sorted(x.unique()))),
            'Contig_Idx': ('Contig_Idx', 'first'),
            'Symbol_Contig_Idx': ('Symbol_Contig_Idx', lambda x: '|'.join(set(sorted(x.unique())))),
            'ntargeted_tx': ('Symbol', 'nunique'),
            'nappearances': ('Symbol', 'count'),
            'exon_id': ('exon_id', lambda x: '|'.join(sorted(x.dropna().unique().astype(str)))),
            'region': ('region', lambda x: '|'.join(sorted(x.dropna().unique().astype(str)))),
            'tags': ('tags', lambda x: '|'.join(sorted(set(
                tag for tags_str in x.dropna() for tag in tags_str.split('|') if tag
            )))),
            'expression_content': ('sum_expression_norm', sum_expression),
            'overlaps_expressed_tx': ('is_expressed', lambda x : x.any()),
            'sum_log_median_expr': ('log_median_expr', sum_expression),
            'sum_log_median_expr_norm': ('log_median_expr_norm', sum_expression),
            'ttm_priority': ('ttm_priority', sum_expression),
            'max_pct_cell_lines_expr': ('pce', 'max'),
            'max_n_cell_lines_expr': ('n_expressed', 'max')
        }
    )

    # 2. Per-Gene Normalization

    # We define the "Max TTM" for a gene as the intensity of its most expressed isoforms.
    # Note: We group by 'Gene' and find the max ttm_priority existing in our collapsed set.
    gene_max_ttm = collapse_df.groupby('Gene')['ttm_priority'].transform('max')
    
    # Rationale: If max_ttm is 0 or -1, the gene is unmeasured/unexpressed. 
    # Otherwise, we scale the guide's score relative to the "best" guide for that gene.
    collapse_df['ttm_gene_norm'] = np.where(
        gene_max_ttm > 0,
        (collapse_df['ttm_priority'] / gene_max_ttm).round(4),
        collapse_df['ttm_priority'] # Preserve -1.0 or 0.0
    )

    return collapse_df

def main():

    parser = argparse.ArgumentParser(description='Collapse')
    parser.add_argument('-g', '--guides_file', type=str, help='path to guides file')
    parser.add_argument('-o', '--output', type=str, help='output file')
    args = parser.parse_args()

    full_data = load_data(args.guides_file)
    coll_data = collapse_guides(full_data)

    # find guides with count greater than expected and remove those from collapsed and non-collapsed files
    # extract these guides and remove them from both collapsed and non-collapsed files
    multiple_self_hits = coll_data[coll_data.nappearances > coll_data.ntargeted_tx].copy()
    guides_to_remove = multiple_self_hits['Guide Sequence'].unique()
    full_data_uniq = full_data[~full_data['Guide Sequence'].isin(guides_to_remove)].copy()
    coll_data_uniq = coll_data[~coll_data['Guide Sequence'].isin(guides_to_remove)].copy()

    full_data_uniq.to_csv(args.output, index=False)
    coll_data_uniq.to_csv(args.output.replace('.csv', '.collapsed.csv'), index=False)

    print(f'[Full]  {full_data.shape[0]:,} guides / {full_data.Gene.nunique():,} genes / {full_data.Symbol.nunique():,} tx')
    print(f'[Coll] {coll_data.shape[0]:,} guides / {coll_data.Gene.nunique():,} genes / {coll_data.Symbol.nunique():,} tx')
    print(f'[Full no mults] {full_data_uniq.shape[0]:,} guides / {full_data_uniq.Gene.nunique():,} genes / {full_data_uniq.Symbol.nunique():,} tx')
    print(f'[Coll no mults] {coll_data_uniq.shape[0]:,} guides / {coll_data_uniq.Gene.nunique():,} genes / {coll_data_uniq.Symbol.nunique():,} tx')

    has_utr = coll_data_uniq.region.str.contains('UTR', na=False, case=False)
    has_cds_lnc = coll_data_uniq.region.str.contains('CDS|lncRNA', na=False, case=False)
    has_canonical = coll_data_uniq.tags.str.contains('canonical', na=False, case=False)

    # UTR-only output: Hits UTR, does NOT hit CDS or lncRNA
    utr_data = coll_data_uniq[has_utr & ~has_cds_lnc].copy()
    utr_data.to_csv(args.output.replace('.csv', '.collapsed.utr.csv'), index=False)
    
    # CDS output: Hits CDS/lncRNA or UTR. 
    # CRITICALLY: ANY guide that hits a UTR must have the canonical tag to be included, 
    # even if it also hits a CDS or lncRNA.
    cds_data = coll_data_uniq[
        (has_cds_lnc | has_utr) & 
        ~(has_utr & ~has_canonical)
    ].copy()
    cds_data.to_csv(args.output.replace('.csv', '.collapsed.cds.csv'), index=False)
    
if __name__ == '__main__':
    main()
