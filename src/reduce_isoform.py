#!/usr/bin/env python
import pandas as pd
import numpy as np
from bisect import bisect_left, insort
from dataclasses import dataclass
import argparse
from typing import List, Dict, Any, Tuple

@dataclass
class SelectionConfig:
    target_n: int = 5
    score_col: str = 'tiger_score'
    pos_col: str = 'position'

class SelectionEngine:
    @staticmethod
    def is_overlapping(pos: int, sorted_list: List[int], min_dist: int) -> bool:
        idx = bisect_left(sorted_list, pos)
        if idx > 0 and pos - sorted_list[idx-1] < min_dist:
            return True
        if idx < len(sorted_list) and sorted_list[idx] - pos < min_dist:
            return True
        return False

    @classmethod
    def select_guides(cls, df: pd.DataFrame, config: SelectionConfig, current_dist: int) -> pd.DataFrame:
        if df.empty: return df
        # Sort by score descending
        sorted_df = df.sort_values(config.score_col, ascending=False).reset_index(drop=True)
        selected_indices, selected_positions = [], []
        scores, starts = sorted_df[config.score_col].values, sorted_df[config.pos_col].values

        for i in range(len(sorted_df)):
            if len(selected_indices) >= config.target_n: break
            # Ensure no overlap with already selected guides
            if not cls.is_overlapping(starts[i], selected_positions, current_dist):
                selected_indices.append(i)
                insort(selected_positions, starts[i])
        
        return sorted_df.iloc[selected_indices]

def parse_position(row, tx_length_dict=None):
    try:
        tids = str(row['transcript_id_group']).split('|')
        primary_tx = max(tids, key=lambda x: tx_length_dict.get(x, 0)) if tx_length_dict else tids[0]
        mapping_str = str(row.get('Symbol_Contig_Idx', ''))
        for item in mapping_str.split('|'):
            if ':' in item:
                tx, pos = item.rsplit(':', 1)
                if tx == primary_tx: return int(pos)
        return row.get('Contig_Idx', -1)
    except: return -1

def process_gene(gene_id: str, gene_df: pd.DataFrame, config: SelectionConfig, target_n: int) -> pd.DataFrame:
    # Filter for guides that target exactly one transcript
    # We look for rows where 'transcript_id_group' has no pipe '|'
    isoform_specific_df = gene_df[~gene_df['transcript_id_group'].str.contains(r'\|', regex=True)].copy()
    
    if isoform_specific_df.empty:
        return pd.DataFrame()

    results = []
    
    # Group by the single transcript ID
    for tid, group_df in isoform_specific_df.groupby('transcript_id_group'):
        # Strategies: We try to get 5 guides. 
        # We can try different distances if we want to be strict, but for now let's just use a reasonable default 
        # or iterate if needed. The prompt says "select not just 1 per gene but all possible".
        # I'll stick to a single pass with a standard distance (e.g. 12bp or similar, matching reduce_constitutive default strategies maybe?)
        # checking reduce_constitutive: it tries 23, 12, 3.
        # Let's try to be strict first (23) then relax? Or just pick one?
        # The prompt says "we want 5 guides per transcript".
        # Let's try to get 5 with strict distance, if not, relax distance? 
        # Simpler approach first: just use a standard distance like 12bp (CRISPR width roughly).
        
        # Actually, let's just use one reasonable distance for now, say 12.
        # If we need more complexity we can add it.
        
        selected = SelectionEngine.select_guides(group_df, config, current_dist=12)
        
        # If we didn't get enough, strictly speaking we just return what we found.
        if not selected.empty:
            # Add metadata
            selected = selected.assign(
                result_class=f"isoform_specific_12bp",
                ntargeted_tx=1
            )
            results.append(selected)

    if results:
        return pd.concat(results)
    return pd.DataFrame()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--metadata")
    parser.add_argument("--output", required=True)
    parser.add_argument("--target_n", type=int, default=5)
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep=None, engine='python')
    # Rename columns to match internal logic if needed, as done in reduce_constitutive
    df = df.rename(columns={'Gene': 'gene_id', 'Symbol': 'transcript_id_group', 'Guide Score': 'tiger_score', 'Title': 'biotype'})

    tx_length_dict = None
    if args.metadata:
        m_df = pd.read_table(args.metadata, low_memory=False)
        tx_length_dict = dict(zip(m_df['transcript_id'].astype(str), m_df['transcript_length']))
    
    df['gene_id'] = df['gene_id'].astype(str).str.strip()
    df['position'] = df.apply(parse_position, tx_length_dict=tx_length_dict, axis=1)
    
    # Ensure columns exist
    for col in ['expression_content', 'tags']:
        if col not in df.columns: df[col] = 0 if col != 'tags' else ""

    config = SelectionConfig(target_n=args.target_n)
    
    all_selected_dfs = []

    # Process by gene (though strictly we could just group by transcript directly, 
    # grouping by gene first keeps structure similar to reduce_constitutive)
    for gene_id, group in df.groupby('gene_id'):
        res_df = process_gene(gene_id, group, config, args.target_n)
        if not res_df.empty:
            all_selected_dfs.append(res_df)

    if all_selected_dfs:
        result_df = pd.concat(all_selected_dfs)
         # Reorder columns as in reduce_constitutive
        cols_to_keep = ['Guide Sequence', 'tiger_score', 'biotype', 'gene_id', 'result_class', 
                        'region', 'position', 'expression_content', 'ntargeted_tx', 
                        'nappearances', 'exon_id', 'transcript_id_group', 'tags']
        # Only keep columns that actually exist
        cols_to_keep = [c for c in cols_to_keep if c in result_df.columns]
        
        result_df = result_df[cols_to_keep].copy()
        result_df.to_csv(args.output, sep="\t", index=False)
    else:
        # Create empty df with columns if possible, or just empty file logic
        # reduce_constitutive doesn't seemingly handle the empty case explicitly for file creation other than doing nothing?
        # Actually if all_selected_dfs is empty it doesn't write anything? 
        # Let's write header at least.
        pd.DataFrame(columns=['Guide Sequence', 'tiger_score', 'biotype', 'gene_id', 'result_class']).to_csv(args.output, sep="\t", index=False)

