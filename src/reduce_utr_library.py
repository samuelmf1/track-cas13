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
from typing import List, Dict, Any, Tuple, Set
import subprocess
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
    def select_guides(cls, df: pd.DataFrame, config: SelectionConfig, current_dist: int, min_score: float, sassy_checker: SassyChecker = None, gene_id: str = None, existing_positions: List[int] = None) -> Tuple[pd.DataFrame, Dict[str, int]]:
        stats = {'total': len(df), 'low_score': 0, 'overlap': 0, 'off_target': 0}
        if df.empty: return df, stats
        sorted_df = df.sort_values(config.score_col, ascending=False).reset_index(drop=True)
        selected_indices, selected_positions = [], sorted(existing_positions) if existing_positions else []
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
    def __init__(self, target_n: int = 5, sassy_checker: SassyChecker = None, utr_data: Dict[str, Dict[str, int]] = None):
        self.target_n = target_n
        self.sassy_checker = sassy_checker
        self.utr_data = utr_data or {}

    def process_gene(self, gene_id: str, gene_df: pd.DataFrame, config: SelectionConfig, canonical_tx: str = None) -> Tuple[pd.DataFrame, List[Dict]]:
        all_tids_raw = gene_df['transcript_id_group'].dropna().unique()
        group_to_tx = {g: set(str(g).split('|')) for g in all_tids_raw}
        all_tx_in_gene = set().union(*group_to_tx.values())

        # Step 1: Prep UTR Data
        tx_to_total_utr = {}
        for tx in all_tx_in_gene:
            u = self.utr_data.get(tx, {'five_prime_UTR': 0, 'three_prime_UTR': 0, 'generic': 0})
            tx_to_total_utr[tx] = u.get('five_prime_UTR', 0) + u.get('three_prime_UTR', 0) + u.get('generic', 0)

        unique_lengths = sorted(set(tx_to_total_utr.values()))
        
        # Step 2: Define and Evaluate Splits
        best_split = None
        best_split_score = -1

        for i in range(len(unique_lengths)):
            for j in range(i, len(unique_lengths)):
                l_short_max = unique_lengths[i]
                l_long_min = unique_lengths[j]
                
                if l_long_min - l_short_max < 50:
                    continue
                
                u_short = {tx for tx, length in tx_to_total_utr.items() if length <= l_short_max}
                u_long = {tx for tx, length in tx_to_total_utr.items() if length >= l_long_min}
                
                # Guides that target ONLY transcripts in u_short
                pool_short = gene_df[gene_df['transcript_id_group'].apply(lambda g: group_to_tx[g].issubset(u_short))]
                # Guides that target ONLY transcripts in u_long
                pool_long = gene_df[gene_df['transcript_id_group'].apply(lambda g: group_to_tx[g].issubset(u_long))]
                
                if pool_short.empty or pool_long.empty:
                    continue

                # Try to select guides
                current_winner_split = None
                for strat in strategies if 'strategies' in locals() else self._get_strategies():
                    sel_short, stats_short = SelectionEngine.select_guides(pool_short, config, strat['dist'], strat['score'], sassy_checker=self.sassy_checker, gene_id=gene_id)
                    if sel_short.empty: continue
                    
                    sel_long, stats_long = SelectionEngine.select_guides(pool_long, config, strat['dist'], strat['score'], sassy_checker=self.sassy_checker, gene_id=gene_id, existing_positions=sel_short[config.pos_col].tolist())
                    if sel_long.empty: continue
                    
                    # Score this split
                    n_total = len(sel_short) + len(sel_long)
                    # Penalize incomplete sets
                    score = n_total
                    if len(sel_short) < self.target_n: score -= 0.5
                    if len(sel_long) < self.target_n: score -= 0.5
                    
                    # Tie-breaker: Tiger score
                    avg_tiger = (sel_short[config.score_col].mean() + sel_long[config.score_col].mean()) / 2
                    score += avg_tiger * 0.1
                    
                    current_winner_split = {
                        'sel_short': sel_short.assign(result_class=strat['label'], targeted_transcript="Short_UTR_Ensemble"),
                        'sel_long': sel_long.assign(result_class=strat['label'], targeted_transcript="Long_UTR_Ensemble"),
                        'label': strat['label'],
                        'score': score,
                        'n_short': len(sel_short),
                        'n_long': len(sel_long),
                        'short_ids': "|".join(sorted(u_short)),
                        'long_ids': "|".join(sorted(u_long)),
                        'fail_reason': ""
                    }
                    break
                
                if current_winner_split and current_winner_split['score'] > best_split_score:
                    best_split = current_winner_split
                    best_split_score = current_winner_split['score']

        summary_rows = []
        biotype = gene_df['biotype'].iloc[0] if 'biotype' in gene_df.columns else "N/A"

        if best_split:
            final_df = pd.concat([best_split['sel_short'], best_split['sel_long']])
            
            # Create summary rows for the two ensembles
            for side in ['Short', 'Long']:
                prefix = side.lower()
                sel_side = best_split[f'sel_{prefix}']
                u_ensemble = best_split[f'{prefix}_ids'].split('|')
                
                # Get metrics from pool
                ensemble_pool = gene_df[gene_df['transcript_id_group'].apply(lambda g: group_to_tx[g].issubset(set(u_ensemble)))]
                
                guide_expression = ensemble_pool['guide_expression'].max() if 'guide_expression' in ensemble_pool.columns else -1
                isexpr = ensemble_pool['overlaps_expressed_tx'].any() if 'overlaps_expressed_tx' in ensemble_pool.columns else False
                guide_expression_norm = ensemble_pool['guide_expression_norm'].max() if 'guide_expression_norm' in ensemble_pool.columns else -1
                max_pct_cell_lines_expr = ensemble_pool['max_pct_cell_lines_expr'].max() if 'max_pct_cell_lines_expr' in ensemble_pool.columns else -1
                max_n_cell_lines_expr = ensemble_pool['max_n_cell_lines_expr'].max() if 'max_n_cell_lines_expr' in ensemble_pool.columns else -1
                
                overlaps_canonical = any(tx == canonical_tx for tx in u_ensemble) if canonical_tx else False
                is_low_expression = guide_expression >= 0 and guide_expression < 1

                summary_rows.append({
                    'gene_id': gene_id,
                    'transcript_id_targeted': f"{side}_UTR_Ensemble",
                    'biotype': biotype,
                    'result_class': best_split['label'],
                    'guide_expression': guide_expression,
                    'guide_expression_norm': guide_expression_norm,
                    'overlaps_expr': isexpr,
                    'n_targeted_guides': best_split[f'n_{prefix}'],
                    'confidence_flag': 'LOW_EXPRESSION' if is_low_expression else 'HIGH_CONFIDENCE',
                    'fail_reason': "",
                    'overlaps_canonical': overlaps_canonical,
                    'max_pct_cell_lines_expr': max_pct_cell_lines_expr,
                    'max_n_cell_lines_expr': max_n_cell_lines_expr,
                    'ensemble_members': best_split[f'{prefix}_ids']
                })
            return final_df, summary_rows
        else:
            # Failure case
            summary_rows.append({
                'gene_id': gene_id,
                'transcript_id_targeted': "No_Qualified_UTR_Split",
                'biotype': biotype,
                'result_class': "NONE",
                'n_targeted_guides': 0,
                'confidence_flag': 'FAILED',
                'fail_reason': "No valid UTR split ≥50bp with double-sided guides"
            })
            return pd.DataFrame(), summary_rows

    def _get_strategies(self):
        score_map = {0.8: 'Q1', 0.6: 'Q2', 0.4: 'Q3', 0.2: 'Q4', 0.0: 'Q5'}
        strategies = []
        idx = 1
        for score in [0.8, 0.6, 0.4, 0.2, 0.0]:
            q_label = score_map[score]
            for dist in [23, 6]:
                strategies.append({'dist': dist, 'score': score, 'label': f"{idx:02d}_{q_label}_{dist}bp", 'idx': idx})
                idx += 1
        return strategies
def load_utr_data(gtf_path: str) -> Dict[str, Dict[str, int]]:
    utr_data = {}
    print(f"Loading UTR data from GTF: {gtf_path}", file=sys.stderr)
    cds_bounds = {}
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS': continue
            match = re.search(r'transcript_id "([^"]+)"', parts[8])
            if match:
                tx = match.group(1)
                start, end, strand = int(parts[3]), int(parts[4]), parts[6]
                if tx not in cds_bounds: cds_bounds[tx] = [strand, start, end]
                else:
                    cds_bounds[tx][1] = min(cds_bounds[tx][1], start)
                    cds_bounds[tx][2] = max(cds_bounds[tx][2], end)

    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] not in ['UTR', 'five_prime_UTR', 'three_prime_UTR']: continue
            match = re.search(r'transcript_id "([^"]+)"', parts[8])
            if match:
                tx = match.group(1)
                start, end, length = int(parts[3]), int(parts[4]), int(parts[4]) - int(parts[3]) + 1
                if tx not in utr_data: utr_data[tx] = {'five_prime_UTR': 0, 'three_prime_UTR': 0, 'generic': 0}
                if parts[2] == 'five_prime_UTR': utr_data[tx]['five_prime_UTR'] += length
                elif parts[2] == 'three_prime_UTR': utr_data[tx]['three_prime_UTR'] += length
                else:
                    if tx in cds_bounds:
                        strand, cds_min, cds_max = cds_bounds[tx]
                        if strand == '+':
                            if end <= cds_min: utr_data[tx]['five_prime_UTR'] += length
                            elif start >= cds_max: utr_data[tx]['three_prime_UTR'] += length
                        else:
                            if start >= cds_max: utr_data[tx]['five_prime_UTR'] += length
                            elif end <= cds_min: utr_data[tx]['three_prime_UTR'] += length
                    else: utr_data[tx]['generic'] += length
    return utr_data

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
    parser.add_argument("--gtf", help="Path to GTF file for UTR length calculation")
    args = parser.parse_args()

    utr_data = load_utr_data(args.gtf) if args.gtf else {}

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
    for col in ['guide_expression', 'overlaps_expressed_tx', 'tags', 'guide_expression_norm', 'max_pct_cell_lines_expr', 'max_n_cell_lines_expr']:
        if col not in df.columns: 
            if col == 'tags': df[col] = ""
            elif col == 'overlaps_expressed_tx': df[col] = False
            else: df[col] = -1
    
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

    selector = GuideSelector(target_n=args.target_n, sassy_checker=sassy_checker, utr_data=utr_data)
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
            'guide_expression', 'guide_expression_norm', 
            'max_pct_cell_lines_expr', 'max_n_cell_lines_expr',
            'overlaps_expressed_tx', 'exon_id', 'transcript_id_group', 'Symbol_Contig_Idx', 'tags']].copy()
        
        guides_selected = result_df[~result_df['targeted_transcript'].isin(excluded_transcripts)]
        guides_excluded = result_df[result_df['targeted_transcript'].isin(excluded_transcripts)]
        
        guides_selected.to_csv(args.output, sep="\t", index=False)
        guides_excluded.to_csv(args.output.replace('.selected.', '.excluded.'), sep="\t", index=False)
    else:
        pd.DataFrame().to_csv(args.output, sep="\t", index=False)
        pd.DataFrame().to_csv(args.output.replace('.selected.', '.excluded.'), sep="\t", index=False)