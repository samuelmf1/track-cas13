#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path
from typing import List, Dict, Optional, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter
from tqdm import tqdm
import re
import io

logger = logging.getLogger(__name__)


def setup_logging(log_file: str):
    """Setup logging to both console and file in the output directory."""
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ],
        force=True  # Clear existing handlers
    )
    logger.info(f"Logging initialized. Log file: {log_file}")

def parse_args():
    p = argparse.ArgumentParser(description='Analyze transcript expression and merge with summary')
    add = p.add_argument
    
    # Required Arguments
    add('--expr', type=str, required=True, help='Path to expression parquet')
    add('--model', type=str, required=True, help='Path to model metadata CSV')
    add('--transcript-summary', type=str, required=True, help='Path to transcript summary TSV')
    
    # Optional Arguments
    add('--output-dir', type=str, default='.', help='Output directory (default: .)')
    add('--tpm-pc', type=float, default=1.0, help='TPM threshold: protein-coding (default: 1.0)')
    add('--tpm-lnc', type=float, default=0.5, help='TPM threshold: lncRNA (default: 0.5)')
    add('--n-cell-lines', type=int, default=17, help='Min cell lines (default: 17)')
    add('--skip-models', action='store_true', help='Skip model processing')

    return p.parse_args()

def load_transcript_summary(tsv_path: str) -> pd.DataFrame:
    """
    Load the transcript summary TSV file (previously generated from GTF).
    This replaces the load_gtf_mapping function.
    """
    logger.info(f"Loading transcript summary file: {tsv_path}")
    
    try:
        mapping_df = pd.read_csv(tsv_path, sep='\t', low_memory=False)
    except Exception as e:
        logger.error(f"Error reading transcript summary file: {e}")
        raise
    
    logger.info(f"Loaded {len(mapping_df):,} transcripts from summary file")
    
    # Report basic stats
    logger.info(f"Total transcripts loaded: {len(mapping_df):,}")
    logger.info(f"Unique genes: {mapping_df['gene_id'].nunique():,}")
    logger.info(f"Unique transcript_ids: {mapping_df['transcript_id'].nunique():,}")
    
    # Report biotype distribution
    if 'biotype' in mapping_df.columns:
        logger.info(f"Biotype distribution:\n{mapping_df['biotype'].value_counts().head(10)}")
    
    # Ensure required columns exist
    required_cols = ['gene_id', 'transcript_id', 'biotype']
    missing_cols = [col for col in required_cols if col not in mapping_df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Rename biotype to gene_type for compatibility with downstream functions
    # (the original script used 'gene_type' internally)
    if 'gene_type' not in mapping_df.columns and 'biotype' in mapping_df.columns:
        mapping_df['gene_type'] = mapping_df['biotype']
    
    return mapping_df

def extract_gencode_version(filepath: str) -> str:
    """
    Extract GENCODE version from filename.
    Looks for patterns like 'gencode.v49', 'v49', 'gencode_v38', etc.
    Returns version string like 'v49' or 'v38' (default if not found).
    """
    filename = Path(filepath).name
    
    # Try various patterns
    patterns = [
        r'gencode[._]?(v\d+)',  # gencode.v49 or gencode_v49
        r'(v\d+)',              # just v49
    ]
    
    for pattern in patterns:
        match = re.search(pattern, filename, re.IGNORECASE)
        if match:
            return match.group(1).lower()
    
    return 'v38'  # Default if no version found

def extract_stable_id(tid: str) -> str:
    """
    Extract stable ID preserving _PAR_Y suffix.
    e.g., ENST00000452144.6 -> ENST00000452144
    e.g., ENST00000452144.6_PAR_Y -> ENST00000452144_PAR_Y
    """
    if '_PAR_Y' in tid:
        # Split off _PAR_Y, remove version, add _PAR_Y back
        base = tid.replace('_PAR_Y', '')
        stable = base.split('.')[0]
        return stable + '_PAR_Y'
    else:
        return tid.split('.')[0]


def analyze_transcript_overlap(omic_df: pd.DataFrame, mapping_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Analyze overlap between expression data and transcript summary annotation.
    Uses stable transcript IDs (without version) for matching since expression data
    may be from a different GENCODE version.
    
    Returns dict with overlap statistics for use in downstream processing.
    """
    # Get transcript IDs from expression data
    omic_txs_full = set([c for c in omic_df.columns if 'ENST' in c])
    omic_txs_stable = set([extract_stable_id(c) for c in omic_txs_full])
    
    # Get transcript IDs from mapping
    mapping_txs_full = set(mapping_df['transcript_id'])
    mapping_txs_stable = set(mapping_df['transcript_id'].apply(extract_stable_id))
    
    # Calculate overlaps using stable IDs
    in_both_stable = omic_txs_stable & mapping_txs_stable
    in_omic_not_summary_stable = omic_txs_stable - mapping_txs_stable
    in_summary_not_omic_stable = mapping_txs_stable - omic_txs_stable
    
    # Also check exact version matches
    in_both_exact = omic_txs_full & mapping_txs_full

    # Find genes that have NO transcripts in the expression data
    # 1. Map all annotation transcripts to stable IDs and check presence in omic data
    temp_df = mapping_df[['gene_id', 'biotype', 'transcript_id']].copy()
    temp_df['stable_id'] = temp_df['transcript_id'].apply(extract_stable_id)
    temp_df['has_match'] = temp_df['stable_id'].isin(omic_txs_stable)
    
    # 2. Group by gene and check if ANY transcript matches
    gene_has_match = temp_df.groupby('gene_id')['has_match'].any()
    
    # 3. Get gene biotypes
    # (Assuming one biotype per gene, take the first)
    gene_biotypes = temp_df.drop_duplicates('gene_id').set_index('gene_id')['biotype']
    
    # 4. Identify genes with NO matches
    missing_genes = gene_has_match[~gene_has_match].index
    
    # 5. Count by biotype
    missing_counts = gene_biotypes.loc[missing_genes].value_counts()
    
    missing_pc_genes = missing_counts.get('protein_coding', 0)
    missing_lnc_genes = missing_counts.get('lncRNA', 0)

    total_genes = mapping_df['gene_id'].nunique()
    protein_coding_genes = mapping_df[mapping_df['biotype'] == 'protein_coding']['gene_id'].nunique()
    lncrna_genes = mapping_df[mapping_df['biotype'] == 'lncRNA']['gene_id'].nunique()
    
    total_transcripts = mapping_df['transcript_id'].nunique()
    protein_coding_txs = mapping_df[mapping_df['biotype'] == 'protein_coding']['transcript_id'].nunique()
    lncrna_txs = mapping_df[mapping_df['biotype'] == 'lncRNA']['transcript_id'].nunique()
    
    # Log analysis
    logger.info("=" * 60)
    logger.info("TRANSCRIPT OVERLAP ANALYSIS")
    logger.info("=" * 60)
    logger.info(f"Expression data transcripts: {len(omic_txs_full):,}")
    logger.info(f"Summary annotation transcripts: {len(mapping_txs_full):,}")
    logger.info("")
    logger.info("Matching on STABLE IDs (without version, e.g., ENST00000456789):")
    logger.info(f"  Transcripts in both: {len(in_both_stable):,}")
    logger.info(f"  In expression only (not in summary): {len(in_omic_not_summary_stable):,}")
    logger.info(f"  In summary only (not in expression): {len(in_summary_not_omic_stable):,}")
    metric = 100 * len(in_both_stable) / len(mapping_txs_stable) if len(mapping_txs_stable) > 0 else 0
    logger.info(f"  Overlap rate (both/summary): {metric:.1f}%")
    logger.info("")
    logger.info("Matching on EXACT IDs (with version, e.g., ENST00000456789.5):")
    logger.info(f"  Exact matches: {len(in_both_exact):,}")
    logger.info(f"  Version mismatches (same locus, different version): {len(in_both_stable) - len(in_both_exact):,}")
    logger.info("")
    logger.info("GENE OVERLAP ANALYSIS (Genes with NO matched transcripts):")
    logger.info(f"  Protein-coding genes with NO transcripts in expression data: {missing_pc_genes:,} / {protein_coding_genes:,} ({100*missing_pc_genes/protein_coding_genes if protein_coding_genes else 0:.1f}%)")
    logger.info(f"  lncRNA genes with NO transcripts in expression data: {missing_lnc_genes:,} / {lncrna_genes:,} ({100*missing_lnc_genes/lncrna_genes if lncrna_genes else 0:.1f}%)")
    logger.info("")
    logger.info("=" * 60)
    logger.info("TRANSCRIPT SUMMARY ANNOTATION")
    logger.info("=" * 60)
    logger.info(f"Total genes: {total_genes:,}")
    logger.info(f"  Protein-coding genes: {protein_coding_genes:,}")
    logger.info(f"  lncRNA genes: {lncrna_genes:,}")
    logger.info(f"Total transcripts: {total_transcripts:,}")
    logger.info(f"  Protein-coding transcripts: {protein_coding_txs:,}")
    logger.info(f"  lncRNA transcripts: {lncrna_txs:,}")
    logger.info("=" * 60)
    
    return {
        'omic_txs_full': omic_txs_full,
        'omic_txs_stable': omic_txs_stable,
        'mapping_txs_stable': mapping_txs_stable,
        'in_both_stable': in_both_stable,
        'in_both_exact': in_both_exact,
    }

def match_by_stable_id(mapping_df: pd.DataFrame, expressed_df: pd.DataFrame) -> pd.DataFrame:
    """
    Handles matching between expression data and annotation using stable IDs.
    Deduplicates based on stable IDs and flags version mismatches.
    """
    mapping_subset = mapping_df[['transcript_id', 'gene_id', 'gene_name', 'biotype', 
                                'transcript_support_level', 'tag']].copy()
    mapping_subset['transcript_id_stable'] = mapping_subset['transcript_id'].apply(extract_stable_id)
    
    # Deduplicate mapping (keeping first version encountered)
    if mapping_subset['transcript_id_stable'].duplicated().any():
        mapping_subset = mapping_subset.drop_duplicates(subset=['transcript_id_stable'], keep='first')

    # Deduplicate expression data
    if expressed_df['transcript_id_stable'].duplicated().any():
        expressed_df = expressed_df.drop_duplicates(subset=['transcript_id_stable'], keep='first')

    v49_lookup = mapping_subset.set_index('transcript_id_stable')['transcript_id']
    
    merged_df = expressed_df.merge(
        mapping_subset.drop(columns=['transcript_id']),
        on='transcript_id_stable',
        how='inner'
    )
    
    merged_df['transcript_id'] = merged_df['transcript_id_stable'].map(v49_lookup)
    
    # Version Mismatch Check
    merged_df['expr_version'] = merged_df['transcript_id_expr'].str.split('.').str[1]
    merged_df['annot_version'] = merged_df['transcript_id'].str.split('.').str[1]
    merged_df['version_mismatch'] = merged_df['expr_version'] != merged_df['annot_version']
    
    return merged_df

import pandas as pd
import numpy as np

def expression(
    mapping_df: pd.DataFrame,
    omic_df: pd.DataFrame,
    biotype: str,
    tpm_threshold: float,
    n_cell_lines: int,
    skip_models: bool = False
) -> pd.DataFrame:
    """
    Filter transcripts by expression and cell line count using a Top-Decile 
    Operational TPM metric to avoid averaging artifacts.
    """
    model_ids = omic_df['ModelID'].copy() if 'ModelID' in omic_df.columns else None
    
    drop_cols = ['Unnamed: 0', 'SequencingID', 'ModelID', 'IsDefaultEntryForModel',
                'ModelConditionID', 'IsDefaultEntryForMC']
    o_df = omic_df.drop(columns=[c for c in drop_cols if c in omic_df.columns])
    
    # DepMap data is usually log2(TPM + 1). Convert back to linear for additive operations.
    vals = o_df.values
    linear_vals = np.power(2, vals) - 1
    
    # 1. Rationale: Sum of linear TPM is additive. 
    # This represents the total "transcript mass" for a gene across all DepMap samples.
    sum_expr = np.sum(linear_vals, axis=0)
    
    # 2. Rationale: Median of the Top 10% (Top Decile).
    # This addresses the "1 high vs 500 low" problem. We ignore the noise of non-expressing 
    # cell lines and focus on the 'operational' level where Cas13 would actually be used.
    # np.percentile handles the thresholding, and we take the median of values above it.
    thresholds = np.percentile(linear_vals, 90, axis=0)
    # Use nanmedian to calculate the typical expression in the "active" lines only.
    masked_vals = np.where(linear_vals >= thresholds, linear_vals, np.nan)
    median_expr = np.nanmedian(masked_vals, axis=0)
    log_median_expr = np.log2(median_expr + 1)

    # 3. Basal expression indicator: >= 1 TPM in >= 1% of cell lines.
    expressed_mask = linear_vals >= tpm_threshold
    n_expr = expressed_mask.sum(axis=0)
    
    # Calculate Percentage of Cell Lines Expressed (PCE)
    pce = (n_expr / linear_vals.shape[0]) * 100
    
    expressed_df = pd.DataFrame({
        "transcript_id_expr": o_df.columns,
        "transcript_id_stable": [extract_stable_id(tid) for tid in o_df.columns],
        "n_expressed": n_expr.astype(int),
        "sum_expression": sum_expr,
        "log_median_expr": log_median_expr,
        "is_expressed": n_expr >= n_cell_lines,
        "pce": pce
    })
    
    matched_df = match_by_stable_id(mapping_df, expressed_df)
    type_df = matched_df[matched_df["biotype"] == biotype].copy()
    
    # 4. Rationale: Ranking by 'sum_expression' (linear).
    # This identifies the dominant isoforms across the entire DepMap dataset.
    type_df["expression_rank"] = (
        type_df.groupby("gene_id")["sum_expression"]
        .rank(method="dense", ascending=False)
        .fillna(999) # Handle non-expressed transcripts
        .astype(int)
    )

    # 5. Rationale: Fraction of Gene Mass (sum_expression_norm).
    # This is your core metric for Cas13. It tells you: "If I target this transcript,
    # what % of the total physical mRNA of this gene am I intercepting?"
    gene_totals_sum = type_df.groupby("gene_id")["sum_expression"].transform('sum')
    type_df["sum_expression_norm"] = np.where(
        gene_totals_sum > 0, type_df["sum_expression"] / gene_totals_sum, 0.0
    )

    # 6. Rationale: Local Relative Abundance (log_median_expr_norm).
    # In the high-expressing lines, what is the relative importance of this isoform?
    # This helps distinguish if an isoform is only dominant in low-expressing lines.
    gene_totals_log_median = type_df.groupby("gene_id")["log_median_expr"].transform('sum')
    type_df["log_median_expr_norm"] = np.where(
        gene_totals_log_median > 0, type_df["log_median_expr"] / gene_totals_log_median, 0.0
    )
    
    # 7. Rationale: Targeted TPM Mass (TTM)
    # This is the "Salient Metric." It combines coverage (efficiency) with 
    # operational intensity (feasibility). 
    # High Score = We target a large % of a highly expressed gene.
    # Low Score = Either we miss the main isoforms, or the gene is barely expressed.
    
    # We use the 'log_median_expr' (Top 10% median) as the weight.
    type_df["ttm_score"] = type_df["sum_expression_norm"] * type_df["log_median_expr"]

    # Optional: Normalize TTM score to a 0-1 scale relative to the whole library 
    # for easier interpretation in reports.
    max_ttm = type_df["ttm_score"].max()
    if max_ttm > 0:
        type_df["ttm_priority"] = type_df["ttm_score"] / max_ttm
    else:
        type_df["ttm_priority"] = 0.0

    # Reset rank for non-expressed transcripts to signify they aren't viable targets.
    type_df.loc[~type_df["is_expressed"], "expression_rank"] = 0
    
    # Model ID mapping for downstream verification of specific cell line hits.
    if not skip_models and model_ids is not None:
        model_id_mapping = {}
        transcript_to_idx = {tid: idx for idx, tid in enumerate(o_df.columns)}
        
        for _, row in type_df.iterrows():
            expr_tid = row['transcript_id_expr']
            if expr_tid in transcript_to_idx:
                col_idx = transcript_to_idx[expr_tid]
                expressing_models = model_ids[expressed_mask[:, col_idx]]
                model_id_mapping[row['transcript_id']] = '|'.join(sorted(set(expressing_models.astype(str))))
        
        type_df["expressing_models"] = type_df["transcript_id"].map(model_id_mapping)
    
    return type_df

def process_expression_results(
        mapping_df: pd.DataFrame, 
        expr_lookup: pd.DataFrame
    ) -> pd.DataFrame:
    """
    Cleans up the merged dataframe, fills missing values, and ensures 
    correct data types for the final output.
    """
    df = mapping_df.copy()
    df = df.set_index('transcript_id').join(expr_lookup, how='left').reset_index()
    
    # Integer and Float fills
    fill_map = {
        'n_expressed': -1,
        'sum_expression': -1.0,
        'expression_rank': 0,
        'sum_expression_norm': -1.0,
        'log_median_expr': -1.0,
        'log_median_expr_norm': -1.0,
        'ttm_priority': -1.0,
        'is_expressed': False,
        'pce': -1,
    }
    for col, val in fill_map.items():
        if col in df.columns:
            df[col] = df[col].fillna(val)

    # Convert specific columns to Nullable Booleans
    # This keeps 'NaN' for transcripts not in the expression dataset
    bool_cols = ['version_mismatch']
    for col in bool_cols:
        if col in df.columns:
            df[col] = df[col].astype('boolean')

    if 'expressing_models' in df.columns:
        df['expressing_models'] = df['expressing_models'].fillna('')

    df['sum_expression'] = round(df['sum_expression'], 6)
    df['sum_expression_norm'] = round(df['sum_expression_norm'], 6)
    df['log_median_expr'] = round(df['log_median_expr'], 6)
    df['log_median_expr_norm'] = round(df['log_median_expr_norm'], 6)
    df['ttm_priority'] = round(df['ttm_priority'], 6)
    df['pce'] = round(df['pce'], 6)

    return df

def main():
    """Main analysis workflow."""
    args = parse_args()
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    input_stem = Path(args.transcript_summary).stem.replace('.summary', '')
    output_log_filename = Path(args.output_dir) / f"{input_stem}.expr.log"
    setup_logging(output_log_filename)

    # Extract GENCODE versions from filenames
    annot_version = extract_gencode_version(args.transcript_summary)
    expr_version = extract_gencode_version(args.expr) if 'gencode' in args.expr.lower() else 'v38'
    
    logger.info("GENCODE versions detected:")
    logger.info(f"  Annotation (transcript summary): {annot_version}")
    logger.info(f"  Expression data: {expr_version} (default if not detected)")
    
    logger.info("Loading expression data...")
    omic_df = pd.read_parquet(args.expr)
    
    logger.info("Loading model metadata...")
    omic_meta_df = pd.read_csv(args.model)
    
    logger.info("Loading transcript summary...")
    mapping_df = load_transcript_summary(args.transcript_summary)
    
    overlap_stats = analyze_transcript_overlap(omic_df, mapping_df)
    
    pc = expression(mapping_df, omic_df, 'protein_coding', args.tpm_pc, args.n_cell_lines, args.skip_models)
    lc = expression(mapping_df, omic_df, 'lncRNA', args.tpm_lnc, args.n_cell_lines, args.skip_models)
    
    expression_df = pd.concat([pc, lc], ignore_index=True)
    
    # Prepare lookup
    expression_cols = ['n_expressed', 'expression_rank', 'sum_expression', 'sum_expression_norm', 'log_median_expr', 'log_median_expr_norm', 'ttm_priority', 'version_mismatch', 'is_expressed', 'pce']
    if not args.skip_models:
        expression_cols.append('expressing_models')
    expr_lookup = expression_df.set_index('transcript_id')[expression_cols]
    
    # Process and Save
    output_df = process_expression_results(mapping_df, expr_lookup)
    
    output_filename = Path(args.output_dir) / f"{input_stem}.expr.tsv"
    output_df.to_csv(output_filename, sep='\t', index=False)
    logger.info(f"Saved to {output_filename}")

if __name__ == "__main__":
    main()