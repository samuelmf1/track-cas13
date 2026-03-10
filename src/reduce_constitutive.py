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
    def select_guides(cls, df: pd.DataFrame, config: SelectionConfig, current_dist: int, min_score: float, max_utr: int, sassy_checker: SassyChecker = None, gene_id: str = None) -> Tuple[pd.DataFrame, Dict[str, int]]:
        stats = {'total': len(df), 'low_score': 0, 'overlap': 0, 'off_target': 0, 'max_utr': 0}
        if df.empty: return df, stats
        sorted_df = df.sort_values(config.score_col, ascending=False).reset_index(drop=True)
        selected_indices, selected_positions = [], []
        scores, starts = sorted_df[config.score_col].values, sorted_df[config.pos_col].values
        sequences = sorted_df['Target Sequence'].values if sassy_checker else None
        regions = sorted_df['region'].values if 'region' in sorted_df.columns else None
        utr_count = 0

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
            
            # Constraint: UTR limit
            is_utr = regions is not None and 'UTR' in str(regions[i])
            if is_utr and utr_count >= max_utr:
                stats['max_utr'] += 1
                continue

            selected_indices.append(i)
            insort(selected_positions, starts[i])
            if is_utr:
                utr_count += 1
        return sorted_df.iloc[selected_indices], stats

class GuideSelector:
    def __init__(self, target_n: int = 5, sassy_checker: SassyChecker = None):
        self.target_n = target_n
        self.sassy_checker = sassy_checker

    def process_gene(self, gene_id: str, gene_df: pd.DataFrame, config: SelectionConfig, canonical_tx: str = None) -> Tuple[pd.DataFrame, Dict]:
        all_groups = {tid: set(str(tid).split('|')) for tid in gene_df['transcript_id_group'].unique()}
        
        # Add Canonical Fallback Group if it doesn't exist as an exact match
        if canonical_tx:
            # Find all guides that target the canonical transcript
            # We filter for rows where 'transcript_id_group' contains the canonical_tx
            # The 'tid' for this group will be a special identifier
            canonical_mask = gene_df['transcript_id_group'].astype(str).str.contains(re.escape(canonical_tx))
            if canonical_mask.any():
                can_fallback_id = f"FALLBACK_{canonical_tx}"
                if can_fallback_id not in all_groups:
                    all_groups[can_fallback_id] = {canonical_tx} # Mark as targeting canonical
        
        all_transcripts = set().union(*[s for s in all_groups.values() if not isinstance(s, str)]) 
        # Note: can_fallback_id is in all_groups keys, its value is {canonical_tx}
        
        biotype = gene_df['biotype'].iloc[0] if 'biotype' in gene_df.columns else "N/A"
        
        score_map = {0.8: 'Q1', 0.6: 'Q2', 0.4: 'Q3', 0.2: 'Q4', 0.0: 'Q5'}
        strategies = []
        idx = 1

        for utr_limit in range(3): # 0, 1, 2
            for score in [0.8, 0.6]:
                q_label = score_map[score]
                for dist in [23, 6]:
                    strategies.append({
                        'max_utr': utr_limit, 'dist': dist, 'score': score, 
                        'label': f"{idx:02d}_{q_label}_{utr_limit}UTR_{dist}bp",
                        'idx': idx
                    })
                    idx += 1
        
        for score in [0.4, 0.2]:
            q_label = score_map[score]
            for utr_limit in range(6):
                for dist in [23, 6]:
                    strategies.append({
                        'max_utr': utr_limit, 'dist': dist, 'score': score, 
                        'label': f"{idx:02d}_{q_label}_{utr_limit}UTR_{dist}bp",
                        'idx': idx
                    })
                    idx += 1
        
        for score in [0.0]:
            q_label = score_map[score]
            for utr_limit in [5]:
                for dist in [3]:
                    strategies.append({
                        'max_utr': utr_limit, 'dist': dist, 'score': score, 
                        'label': f"{idx:02d}_{q_label}_{utr_limit}UTR_{dist}bp",
                        'idx': idx
                    })
                    idx += 1

        best_idx_found = None
        candidates = []
        best_failure_stats = None

        # Iterate through strategies to find the first valid rank and look ahead
        for strat in strategies:
            if best_idx_found is not None and strat['idx'] > best_idx_found + config.lookahead:
                break
            
            current_strat_best_stats = None

            for tid, primary_set in all_groups.items():
                if str(tid).startswith("FALLBACK_"):
                    # For fallback, we take all guides that target the canonical transcript
                    orig_can_tx = tid.replace("FALLBACK_", "")
                    group_df = gene_df[gene_df['transcript_id_group'].astype(str).str.contains(re.escape(orig_can_tx))]
                else:
                    group_df = gene_df[gene_df['transcript_id_group'] == tid]
                
                # Prefer CDS/lncRNA/UTR guides, but fall back to all guides if insufficient
                pool = group_df[group_df['region'].astype(str).str.contains('CDS|lncRNA|UTR', na=False)]
                if len(pool) < config.target_n:
                    pool = group_df
                
                selected, stats = SelectionEngine.select_guides(pool, config, strat['dist'], strat['score'], strat['max_utr'], self.sassy_checker, gene_id)

                if current_strat_best_stats is None or stats['total'] > current_strat_best_stats['total']:
                    current_strat_best_stats = stats
                
                if len(selected) == self.target_n:
                    if best_idx_found is None:
                        best_idx_found = strat['idx']
                    
                    # Calculate target set for summary reporting
                    target_set = primary_set
                    target_group_str = "|".join(sorted(target_set))

                    guide_expression = group_df['guide_expression'].iloc[0] if 'guide_expression' in group_df.columns else -1
                    isexpr = group_df['overlaps_expressed_tx'].iloc[0] if 'overlaps_expressed_tx' in group_df.columns else False
                    guide_expression_norm = group_df['guide_expression_norm'].iloc[0] if 'guide_expression_norm' in group_df.columns else -1
                    pce = group_df['max_pct_cell_lines_expr'].iloc[0] if 'max_pct_cell_lines_expr' in group_df.columns else -1
                    n_expressed = group_df['max_n_cell_lines_expr'].iloc[0] if 'max_n_cell_lines_expr' in group_df.columns else -1
                    
                    is_fallback = str(tid).startswith("FALLBACK_")
                    if is_fallback:
                        tags = "Ensembl_canonical_fallback"
                    else:
                        tags = group_df['tags'].iloc[0] if 'tags' in group_df.columns and pd.notna(group_df['tags'].iloc[0]) else "None"
                    overlaps_canonical = canonical_tx in primary_set if canonical_tx else ('canonical' in tags.lower())

                    candidates.append({
                        'df': selected.assign(result_class=strat['label']),
                        'guide_expression': guide_expression, 
                        'guide_expression_norm': guide_expression_norm, 
                        'overlaps_expr':isexpr, 'n_tx': len(primary_set), 'rank': strat['idx'],
                        'tid': tid, 'label': strat['label'], 
                        'biotype': group_df['biotype'].iloc[0] if 'biotype' in group_df.columns else "N/A",
                        'target_group_str': "Canonical_Fallback" if is_fallback else target_group_str, 
                        'target_set': primary_set,
                        'tags': tags, 'overlaps_canonical': overlaps_canonical,
                        'is_fallback': is_fallback,
                        'fail_reason': None,
                        'pce': pce,
                        'n_expressed': n_expressed
                    })
            
            if current_strat_best_stats:
                best_failure_stats = current_strat_best_stats

        if candidates:
            # Sorting logic:
            # 1. Use guide_expression primarily.
            # 2. UNLESS expression is low (norm < 1), then prioritize canonical.
            gene_is_high_expr = any(c['guide_expression'] >= 1 for c in candidates)
            
            if gene_is_high_expr:
                candidates.sort(key=lambda x: x['guide_expression'], reverse=True)
            else:
                candidates.sort(key=lambda x: (x['overlaps_canonical'], x['guide_expression']), reverse=True)
            
            winner = candidates[0]
            
            # Aggregate Symbol_Contig_Idx for the selected guides
            # These are the mappings for each guide across all targeted transcripts
            mapping_list = []
            if 'Symbol_Contig_Idx' in winner['df'].columns:
                mapping_list = winner['df']['Symbol_Contig_Idx'].dropna().astype(str).tolist()
            mapping_summary = ";".join(mapping_list) if mapping_list else "None"

            is_low_expression = winner['guide_expression'] < 1 # 1 TPM

            # If fallback was used, the targeted group reported should ideally be the canonical transcript only
            # but we'll stick to the winner['target_group_str'] which we set above.

            summary_info = {
                'gene_id': gene_id, 'biotype': winner['biotype'], 
                'transcript_id_group_targeted': winner['target_group_str'],
                'transcript_id_mapping_targeted': mapping_summary,
                'transcript_ids_not_targeted': "|".join(all_transcripts - winner['target_set']),
                'result_class': winner['label'], 'guide_expression': winner['guide_expression'], 
                'guide_expression_norm': winner['guide_expression_norm'], 'overlaps_expressed_tx': winner['overlaps_expr'], 'n_targeted_tx': len(winner['target_set']),
                'tags': winner['tags'],
                'confidence_flag': 'LOW_EXPRESSION' if winner['guide_expression'] < 1 else 'HIGH_CONFIDENCE',
                'fail_reason': None,
                'overlaps_canonical': winner['overlaps_canonical'],
                'max_pct_cell_lines_expr': winner['pce'],
                'max_n_cell_lines_expr': winner['n_expressed']
            }
            return winner['df'], summary_info

        # Determine fail reason from stats
        fail_reason = f"< {self.target_n} guides"
        if best_failure_stats:
            if best_failure_stats['total'] < self.target_n:
                fail_reason = f"Insufficient candidates ({best_failure_stats['total']})"
            else:
                reasons = {k: v for k, v in best_failure_stats.items() if k != 'total'}
                if reasons:
                    max_reason = max(reasons, key=reasons.get)
                    count = reasons[max_reason]
                    human_reasons = {'low_score': 'Low Score', 'overlap': 'Overlap', 'off_target': 'Off-target', 'max_utr': 'Max UTR'}
                    fail_reason = f"{human_reasons.get(max_reason, max_reason)} ({count})"

        biotype = gene_df['biotype'].iloc[0] if 'biotype' in gene_df.columns else "N/A"
        return pd.DataFrame(), {
            'gene_id': gene_id, 'biotype': biotype, 
            'transcript_id_group_targeted': "None",
            'transcript_id_mapping_targeted': "None",
            'transcript_ids_not_targeted': "|".join(sorted(all_transcripts)),
            'result_class': "999_failed", 'guide_expression': -1, 'guide_expression_norm': -1, 'overlaps_expressed_tx': -1, 'n_targeted_tx': 0, 'tags': "None",
            'confidence_flag': 'LOW_CONFIDENCE',
            'fail_reason': fail_reason,
            'overlaps_canonical': False,
            'max_pct_cell_lines_expr': -1,
            'max_n_cell_lines_expr': -1
        }

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
    
    # Ensure columns exist (mapping for old and new)
    for col in ['guide_expression', 'overlaps_expressed_tx', 'tags', 'guide_expression_norm', 'expression_rank', 'version_mismatch', 'max_pct_cell_lines_expr', 'max_n_cell_lines_expr']:
        if col not in df.columns: 
            if col == 'version_mismatch': df[col] = False
            elif col == 'overlaps_expressed_tx': df[col] = False
            elif col == 'tags': df[col] = ""
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

    selector = GuideSelector(target_n=args.target_n, sassy_checker=sassy_checker)
    all_selected_dfs, summary_rows = [], []

    for gene_id, group in df.groupby('gene_id'):
        can_tx = canonical_map.get(gene_id)
        res_df, summary_dict = selector.process_gene(gene_id, group, config, canonical_tx=can_tx)
        summary_rows.append(summary_dict)
        if not res_df.empty:
            all_selected_dfs.append(res_df)

    summary_df = pd.DataFrame(summary_rows)

    # Apply exclusion reason logic
    # Only exclude if:
    #   1. Gene failed entirely (999_failed), OR
    #   2. Gene has a known canonical transcript AND the winner doesn't overlap it
    # Genes with no canonical metadata are KEPT.
    def get_exclusion_reason(row):
        if row['result_class'] == '999_failed':
            return row.get('fail_reason', 'No guides')
        if row['gene_id'] in canonical_map and not row.get('overlaps_canonical', False):
            return 'Non-canonical'
        return ''

    summary_df['exclusion_reason'] = summary_df.apply(get_exclusion_reason, axis=1)

    # Filtering logic for excluded genes
    excluded_mask = summary_df['exclusion_reason'] != ''
    excluded_genes = set(summary_df[excluded_mask]['gene_id'].unique())

    # Split and write summaries
    summary_selected = summary_df[~summary_df['gene_id'].isin(excluded_genes)]
    summary_excluded = summary_df[summary_df['gene_id'].isin(excluded_genes)]
    
    # Drop fail_reason helper column if desired, or keep it. keeping it is harmless but we might want to clean up.
    # We'll explicitly select columns or just write all.
    
    summary_selected.to_csv(args.output.replace('.tsv', '') + '.summary.tsv', sep="\t", index=False)
    summary_excluded.to_csv(args.output.replace('.selected.', '.excluded.').replace('.tsv', '') + '.summary.tsv', sep="\t", index=False)

    # Split and write guides
    if all_selected_dfs:
        result_df = pd.concat(all_selected_dfs)
        result_df = result_df[
            ['Guide Sequence', 'tiger_score', 'biotype', 'gene_id', 'result_class', 'region', 
            'guide_expression', 'guide_expression_norm',
            'max_pct_cell_lines_expr', 'max_n_cell_lines_expr', 'expression_rank', 'version_mismatch',
            'overlaps_expressed_tx', 'exon_id', 'transcript_id_group', 'Symbol_Contig_Idx', 'tags']].copy()
        
        guides_selected = result_df[~result_df['gene_id'].isin(excluded_genes)]
        guides_excluded = result_df[result_df['gene_id'].isin(excluded_genes)]
        
        guides_selected.to_csv(args.output, sep="\t", index=False)
        guides_excluded.to_csv(args.output.replace('.selected.', '.excluded.'), sep="\t", index=False)
    else:
        # Create empty Files if no guides at all
        pd.DataFrame().to_csv(args.output, sep="\t", index=False)
        pd.DataFrame().to_csv(args.output.replace('.selected.', '.excluded.'), sep="\t", index=False)