#!/usr/bin/env python
import pandas as pd
import numpy as np
from bisect import bisect_left, bisect_right, insort
from dataclasses import dataclass
import argparse
import os
import sys
import time
import re
from typing import List, Dict, Any, Tuple
try:
    import sassy
    import pysam
except ImportError:
    sassy = None
    pysam = None

@dataclass
class SelectionConfig:
    target_n: int = 5
    score_col: str = 'tiger_score'
    pos_col: str = 'position'
    lookahead: int = 0
    sassy_ref: str = None
    sassy_cache: str = None

class SassyChecker:
    CIGAR_RE = re.compile(r"(\d*)([=XID])")

    def __init__(self, ref_fasta, cache_file=None):
        if not sassy or not pysam:
            raise ImportError("sassy or pysam not installed. Please install sassy-rs and pysam.")
        self.searcher = sassy.Searcher("iupac")
        self.text, self.intervals, self.starts = self._load_ref(ref_fasta)
        self.cache = set()
        self.cache_file = cache_file
        if cache_file and os.path.exists(cache_file):
            with open(cache_file, 'r') as f:
                for line in f:
                    self.cache.add(line.strip())

    def _load_ref(self, fasta_path):
        intervals = []
        starts = []
        
        # We'll build the text dynamically if needed, or pre-allocate and truncate.
        # But bytearray is efficient for append or pre-allocate.
        with pysam.FastaFile(fasta_path) as fasta:
            text = bytearray()
            current_offset = 0
            for name in fasta.references:
                try:
                    seq_str = fasta.fetch(name)
                    if seq_str is None:
                        print(f"Warning: Could not fetch sequence for {name}. Skipping.", file=sys.stderr)
                        continue
                    seq = seq_str.upper().encode('ascii')
                except (ValueError, KeyError, IndexError, OSError) as e:
                    print(f"Warning: Error fetching sequence for {name}: {e}. Skipping.", file=sys.stderr)
                    continue

                parts = name.split('|')
                gene_id = parts[1] if len(parts) > 1 else name
                
                starts.append(current_offset)
                intervals.append((current_offset, current_offset + len(seq), gene_id))
                
                # Append to bytearray
                text.extend(seq)
                text.append(ord(b'N'))
                current_offset += len(seq) + 1
                
        return bytes(text), intervals, starts

    def is_disqualified(self, pattern, gene_id):
        if pattern in self.cache:
            return True
        
        # Search with k=2 to allow for mismatches/indels
        matches = self.searcher.search(pattern.encode('ascii'), self.text, k=2)
        for m in matches:
            if m.strand != '+':
                continue
            
            # Skip matches that include wildcards (Ns) in the reference
            ref_seq = self.text[m.text_start:m.text_end]
            if b'N' in ref_seq:
                continue

            is_strong = self._is_strong_match(m.cigar)
            if is_strong:
                matched_gene_id = self._get_gene_id(m.text_start)
                if matched_gene_id != gene_id:
                    if self.cache_file:
                        with open(self.cache_file, 'a') as f:
                            f.write(pattern + '\n')
                    self.cache.add(pattern)
                    return True
        return False

    def _is_strong_match(self, cigar: str) -> bool:
        """
        Heuristic check for a 'strong' match based on CIGAR string.
        Derived from filter_cigar.py.
        """
        counts = {"=": 0, "X": 0, "I": 0, "D": 0}
        positions = []
        ref_pos = 1

        for length_str, op in self.CIGAR_RE.findall(cigar):
            length = int(length_str) if length_str else 1
            start = ref_pos
            end = ref_pos + (length - 1 if op in ("=", "X", "D") else 0)
            positions.append((op, start, end))
            if op in ("=", "X", "D"):
                ref_pos += length
            counts[op] += length

        total_ref_length = sum(counts.values()) - counts["I"]
        if total_ref_length != 23:
            return False

        x, i_count, d_count = counts["X"], counts["I"], counts["D"]
        if x + i_count + d_count > 2:
            return False

        # Rule 1: up to 2 mismatches, no indels
        if x <= 2 and i_count == 0 and d_count == 0:
            return True

        # Rule 2: any I/D, no mismatches, must be valid by positions
        if x == 0 and (i_count > 0 or d_count > 0):
            # valid_indel_positions logic
            i = 0
            n = len(positions)
            while i < n:
                op, start, _ = positions[i]
                if op in ("I", "D"):
                    run_count = 1
                    j = i + 1
                    while j < n and positions[j][0] == op:
                        run_count += 1
                        j += 1
                    if run_count > 2: return False
                    if start > 3 and start < 21: return False
                    i = j
                else:
                    i += 1
            return True
        
        return False

    def _get_gene_id(self, offset):
        idx = bisect_right(self.starts, offset) - 1
        if idx >= 0:
            start, end, gid = self.intervals[idx]
            if start <= offset < end:
                return gid
        return None

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
    def select_guides(cls, df: pd.DataFrame, config: SelectionConfig, current_dist: int, min_score: float, sassy_checker: SassyChecker = None, gene_id: str = None) -> Tuple[pd.DataFrame, Dict[str, int]]:
        stats = {'total': len(df), 'low_score': 0, 'overlap': 0, 'off_target': 0}
        if df.empty: return df, stats
        sorted_df = df.sort_values(config.score_col, ascending=False).reset_index(drop=True)
        selected_indices, selected_positions = [], []
        scores, starts = sorted_df[config.score_col].values, sorted_df[config.pos_col].values
        sequences = sorted_df['Target Sequence'].values if sassy_checker else None

        for i in range(len(sorted_df)):
            if len(selected_indices) >= config.target_n: break
            
            if scores[i] < min_score:
                stats['low_score'] += 1
                continue

            if cls.is_overlapping(starts[i], selected_positions, current_dist):
                stats['overlap'] += 1
                continue

            if sassy_checker and sassy_checker.is_disqualified(sequences[i], gene_id):
                stats['off_target'] += 1
                continue
            
            selected_indices.append(i)
            insort(selected_positions, starts[i])
        return sorted_df.iloc[selected_indices], stats

class GuideSelector:
    def __init__(self, target_n: int = 5, sassy_checker: SassyChecker = None):
        self.target_n = target_n
        self.sassy_checker = sassy_checker

    def process_gene(self, gene_id: str, gene_df: pd.DataFrame, config: SelectionConfig, canonical_tx: str = None) -> Tuple[pd.DataFrame, List[Dict]]:
        all_tids_raw = gene_df['transcript_id_group'].dropna().unique()
        all_transcripts = set()
        for t in all_tids_raw:
            all_transcripts.update(str(t).split('|'))
        
        biotype = gene_df['biotype'].iloc[0] if 'biotype' in gene_df.columns else "N/A"
        
        score_map = {0.8: 'Q1', 0.6: 'Q2', 0.4: 'Q3', 0.2: 'Q4', 0.0: 'Q5'}
        strategies = []
        idx = 1

        for score in [0.8, 0.6]:
            q_label = score_map[score]
            for dist in [23, 6]:
                strategies.append({'dist': dist, 'score': score, 'label': f"{idx:02d}_{q_label}_{dist}bp", 'idx': idx})
                idx += 1
        
        for score in [0.4, 0.2, 0.0]:
            q_label = score_map[score]
            for dist in [23, 6]:
                strategies.append({'dist': dist, 'score': score, 'label': f"{idx:02d}_{q_label}_{dist}bp", 'idx': idx})
                idx += 1

        selected_dfs = []
        summary_rows = []

        # We want to target every unique combination of transcripts (the powerset)
        all_groups = gene_df['transcript_id_group'].dropna().unique()

        for group_str in sorted(all_groups):
            # Select guides that target exactly this combination of transcripts
            pool = gene_df[gene_df['transcript_id_group'] == group_str]
            
            winner = None
            best_failure_stats = None

            for strat in strategies:
                selected, stats = SelectionEngine.select_guides(pool, config, strat['dist'], strat['score'], sassy_checker=self.sassy_checker, gene_id=gene_id)
                
                if best_failure_stats is None or stats['total'] > best_failure_stats['total']:
                    best_failure_stats = stats
                
                if len(selected) == self.target_n:
                    winner = {
                        'df': selected.assign(result_class=strat['label'], targeted_transcript=group_str),
                        'n_guides': len(selected),
                        'label': strat['label'],
                        'fail_reason': ""
                    }
                    break
            
            if not winner:
                strat = strategies[-1]
                selected, stats = SelectionEngine.select_guides(pool, config, strat['dist'], strat['score'], sassy_checker=self.sassy_checker, gene_id=gene_id)
                
                fail_reason = f"< {self.target_n} guides"
                if best_failure_stats:
                    if best_failure_stats['total'] < self.target_n:
                        fail_reason = f"Insufficient candidates ({best_failure_stats['total']})"
                    else:
                        reasons = {k: v for k, v in best_failure_stats.items() if k != 'total'}
                        if reasons:
                            max_reason = max(reasons, key=reasons.get)
                            count = reasons[max_reason]
                            human_reasons = {'low_score': 'Low Score', 'overlap': 'Overlap', 'off_target': 'Off-target'}
                            fail_reason = f"{human_reasons.get(max_reason, max_reason)} ({count})"
                
                winner = {
                    'df': selected.assign(result_class=f"{strat['label']}_INCOMPLETE", targeted_transcript=group_str),
                    'n_guides': len(selected),
                    'label': f"{strat['label']}_INCOMPLETE",
                    'fail_reason': fail_reason
                }
            
            if not winner['df'].empty:
                selected_dfs.append(winner['df'])
            
            expr = pool['expression_content'].max() if 'expression_content' in pool.columns else -1
            isexpr = pool['overlaps_expressed_tx'].any() if 'overlaps_expressed_tx' in pool.columns else False
            sum_log_median_expr_norm = pool['sum_log_median_expr_norm'].max() if 'sum_log_median_expr_norm' in pool.columns else -1 
            ttm_priority = pool['ttm_priority'].max() if 'ttm_priority' in pool.columns else -1
            ttm_gene_norm = pool['ttm_gene_norm'].max() if 'ttm_gene_norm' in pool.columns else -1
            max_pct_cell_lines_expr = pool['max_pct_cell_lines_expr'].max() if 'max_pct_cell_lines_expr' in pool.columns else -1
            max_n_cell_lines_expr = pool['max_n_cell_lines_expr'].max() if 'max_n_cell_lines_expr' in pool.columns else -1
            
            is_low_expression = sum_log_median_expr_norm >= 0 and sum_log_median_expr_norm < 1
            
            tx_list = set(group_str.split('|'))
            overlaps_canonical = (canonical_tx in tx_list) if canonical_tx else False
            
            summary_info = {
                'gene_id': gene_id, 
                'transcript_id_targeted': group_str,
                'biotype': biotype, 
                'result_class': winner['label'], 
                'expression': expr, 
                'sum_log_median_expr_norm': sum_log_median_expr_norm, 
                'ttm_priority': ttm_priority,
                'ttm_gene_norm': ttm_gene_norm, 
                'overlaps_expr': isexpr, 
                'n_targeted_guides': winner['n_guides'],
                'confidence_flag': 'LOW_EXPRESSION' if is_low_expression else 'HIGH_CONFIDENCE',
                'fail_reason': winner['fail_reason'],
                'overlaps_canonical': overlaps_canonical,
                'max_pct_cell_lines_expr': max_pct_cell_lines_expr,
                'max_n_cell_lines_expr': max_n_cell_lines_expr
            }
            summary_rows.append(summary_info)

        if selected_dfs:
            return pd.concat(selected_dfs), summary_rows
        
        return pd.DataFrame(), summary_rows

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
    except Exception: return -1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv", required=True)
    parser.add_argument("--metadata")
    parser.add_argument("--output", required=True)
    parser.add_argument("--target_n", type=int, default=5)
    parser.add_argument("--lookahead", type=int, default=0)
    parser.add_argument("--sassy_ref", help="Path to reference FASTA for off-target check")
    parser.add_argument("--sassy_cache", help="Path to cache file for disqualified rows")
    args = parser.parse_args()

    df = pd.read_csv(args.tsv, sep=None, engine='python')
    df = df.rename(columns={'Gene': 'gene_id', 'Symbol': 'transcript_id_group', 'Guide Score': 'tiger_score', 'Title': 'biotype'})

    tx_length_dict = None
    canonical_map = {}
    if args.metadata:
        m_df = pd.read_table(args.metadata, low_memory=False)
        tx_length_dict = dict(zip(m_df['transcript_id'].astype(str), m_df['transcript_length']))
        if 'is_representative' in m_df.columns:
            can_df = m_df[m_df['is_representative'] == True]
            canonical_map = dict(zip(can_df['gene_id'].astype(str), can_df['transcript_id'].astype(str)))
        elif 'tag' in m_df.columns:
            can_df = m_df[m_df['tag'].str.contains('Ensembl_canonical', na=False)]
            canonical_map = dict(zip(can_df['gene_id'].astype(str), can_df['transcript_id'].astype(str)))

    df['gene_id'] = df['gene_id'].astype(str).str.strip()
    df['position'] = df.apply(parse_position, tx_length_dict=tx_length_dict, axis=1)
    for col in ['expression_content', 'overlaps_expressed_tx', 'tags']:
        if col not in df.columns: df[col] = 0 if col != 'tags' else ""
    
    config = SelectionConfig(
        target_n=args.target_n, lookahead=args.lookahead, 
        sassy_ref=args.sassy_ref, sassy_cache=args.sassy_cache
    )
    
    sassy_checker = None
    if args.sassy_ref:
        print(f"Loading Sassy reference: {args.sassy_ref}", file=sys.stderr)
        start_time = time.time()
        sassy_checker = SassyChecker(args.sassy_ref, args.sassy_cache)
        print(f"Loaded Sassy reference in {time.time() - start_time:.2f}s", file=sys.stderr)

    selector = GuideSelector(target_n=args.target_n, sassy_checker=sassy_checker)
    all_selected_dfs, summary_rows = [], []

    for gene_id, group in df.groupby('gene_id'):
        can_tx = canonical_map.get(gene_id)
        res_df, summary_dicts = selector.process_gene(gene_id, group, config, canonical_tx=can_tx)
        summary_rows.extend(summary_dicts)
        if not res_df.empty:
            all_selected_dfs.append(res_df)

    summary_df = pd.DataFrame(summary_rows)

    def get_exclusion_reason(row):
        if row['n_targeted_guides'] == 0:
            return row.get('fail_reason', 'No guides')
        if row.get('fail_reason') and row['fail_reason'] != "":
            return row['fail_reason']
        return ''

    summary_df['exclusion_reason'] = summary_df.apply(get_exclusion_reason, axis=1)

    excluded_mask = summary_df['exclusion_reason'] != ''
    excluded_transcripts = set(summary_df[excluded_mask]['transcript_id_targeted'].unique())

    summary_selected = summary_df[~summary_df['transcript_id_targeted'].isin(excluded_transcripts)]
    summary_excluded = summary_df[summary_df['transcript_id_targeted'].isin(excluded_transcripts)]
    
    summary_selected.to_csv(args.output.replace('.tsv', '') + '.summary.tsv', sep="\t", index=False)
    summary_excluded.to_csv(args.output.replace('.selected.', '.excluded.').replace('.tsv', '') + '.summary.tsv', sep="\t", index=False)

    if all_selected_dfs:
        result_df = pd.concat(all_selected_dfs)
        result_df = result_df[
            ['Guide Sequence', 'tiger_score', 'biotype', 'gene_id', 'result_class', 'region', 'targeted_transcript',
            'expression_content', 'sum_log_median_expr_norm', 'ttm_priority', 'ttm_gene_norm', 
            'max_pct_cell_lines_expr', 'max_n_cell_lines_expr',
            'overlaps_expressed_tx', 'exon_id', 'transcript_id_group', 'Symbol_Contig_Idx', 'tags']].copy()
        
        guides_selected = result_df[~result_df['targeted_transcript'].isin(excluded_transcripts)]
        guides_excluded = result_df[result_df['targeted_transcript'].isin(excluded_transcripts)]
        
        guides_selected.to_csv(args.output, sep="\t", index=False)
        guides_excluded.to_csv(args.output.replace('.selected.', '.excluded.'), sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(args.output, sep="\t", index=False)
        pd.DataFrame().to_csv(args.output.replace('.selected.', '.excluded.'), sep="\t", index=False)