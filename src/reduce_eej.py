#!/usr/bin/env python
import pandas as pd
import numpy as np
from bisect import bisect_left, insort
from dataclasses import dataclass
import argparse
from typing import List, Dict, Any, Tuple
import re

@dataclass
class SelectionConfig:
    target_n: int = 5
    score_col: str = 'tiger_score'
    pos_col: str = 'position'
    min_dist: int = 0  # Default minimum distance if not specified

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
    def select_guides(cls, df: pd.DataFrame, config: SelectionConfig) -> pd.DataFrame:
        if df.empty: return df
        # Sort by score descending
        sorted_df = df.sort_values(config.score_col, ascending=False).reset_index(drop=True)
        selected_indices, selected_positions = [], []
        scores, starts = sorted_df[config.score_col].values, sorted_df[config.pos_col].values

        for i in range(len(sorted_df)):
            if len(selected_indices) >= config.target_n: break
            # Ensure no overlap
            if not cls.is_overlapping(starts[i], selected_positions, config.min_dist):
                selected_indices.append(i)
                insort(selected_positions, starts[i])
        
        return sorted_df.iloc[selected_indices]

def validate_exon_id(exon_id: str) -> bool:
    """
    Validates that exon_id matches format ENSE....-ENSE.... 
    and does NOT contain '_inadmissible'.
    """
    if not isinstance(exon_id, str): return False
    if '_inadmissible' in exon_id: return False
    # Regex for ENSE<digits>.<version>-ENSE<digits>.<version>
    # Note: version is optional in some contexts but user example had .1
    # User said: ENSE00003384789.1-ENSE00003662717.1
    # I'll use a slightly more permissive regex for the ID part to be safe, but strict on structure.
    # ^ENSE\d+(\.\d+)?-ENSE\d+(\.\d+)?$
    pattern = r'^ENSE\d+(\.\d+)?-ENSE\d+(\.\d+)?$'
    return bool(re.match(pattern, exon_id))

def parse_position(row):
    # reduce_constitutive.py uses a complex parse_position. 
    # For EEJ, we might just use 'Contig_Idx' or similar if it exists, or verify what 'position' means.
    # The user didn't specify position logic changes, so I will try to adapt the existing one 
    # but simplify it since we know there is only 1 transcript.
    try:
        mapping_str = str(row.get('Symbol_Contig_Idx', ''))
        # mapping_str format: TRANSCRIPT:POS|...
        # Since we filter for single transcript group, we expect 1 relevant entry.
        # But wait, 'transcript_id_group' has the ID.
        tid = row['transcript_id_group']
        for item in mapping_str.split('|'):
            if ':' in item:
                tx, pos = item.rsplit(':', 1)
                if tx == tid: return int(pos)
        return row.get('Contig_Idx', -1)
    except: return -1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--target_n", type=int, default=5)
    # user didn't extensively specify metadata input or tx lengths, so I'll omit for now unless needed for position.
    # reduce_constitutive used metadata for tx length to pick primary. Here we only have 1 tx per group.
    
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep=None, engine='python')
    
    # Standardize column names if needed (copying from reduce_constitutive)
    rename_map = {
        'Gene': 'gene_id', 
        'Symbol': 'transcript_id_group', 
        'Guide Score': 'tiger_score', 
        'Title': 'biotype'
    }
    df = df.rename(columns={k:v for k,v in rename_map.items() if k in df.columns})

    # Filter 1: Single transcript in group
    # User correction: guides CAN target multiple transcripts if they share the exon-exon junction.
    # So we do NOT filter out groups with pipes.
    # df = df[~df['transcript_id_group'].astype(str).str.contains(r'\|')]

    
    # Filter 2: Exon ID format
    if 'exon_id' in df.columns:
        df = df[df['exon_id'].apply(validate_exon_id)]
    
    # Calculate position
    # The 'parse_position' needs to work. 
    # If 'position' col doesn't exist or we want to re-parse.
    # reduce_constitutive parses it.
    df['position'] = df.apply(parse_position, axis=1) # Simplified version

    # Config
    config = SelectionConfig(target_n=args.target_n)
    
    all_selected_dfs = []
    
    # Strategy: Group by transcript_id_group (which is unique per row now)
    # actually we can just group by transcript_id_group directly.
    # "select not just 1 transcript group per gene but all possible groups."
    
    for tid, group in df.groupby('transcript_id_group'):
        selected = SelectionEngine.select_guides(group, config)
        if not selected.empty:
            all_selected_dfs.append(selected)

    if all_selected_dfs:
        result_df = pd.concat(all_selected_dfs)
        # Ensure columns are present before saving
        cols_to_save = ['Guide Sequence', 'tiger_score', 'biotype', 'gene_id', 'position', 
                        'exon_id', 'transcript_id_group']
        # Add others if they exist
        for c in ['expression_content', 'tags', 'result_class', 'region', 'ntargeted_tx', 'nappearances']:
            if c in df.columns: cols_to_save.append(c)
            
        # Filter cols that actually exist in result_df
        final_cols = [c for c in cols_to_save if c in result_df.columns]
        result_df = result_df[final_cols]
        
        result_df.to_csv(args.output, sep="\t", index=False)
    else:
        # Create empty file with headers if possible or just touch
        # But better to write header of input df if available
        # logic: if input has cols, write empty with those cols
        # if not, nothing.
        with open(args.output, 'w') as f:
            f.write("\t".join(df.columns) + "\n")

