#!/usr/bin/env python3
"""
Validation script for FASTA reformatter output.

Verifies that the output FASTA and summary TSV meet all specifications:
- Header format
- Biotype filtering
- Sequence validity
- Length requirements
- No duplicates
- FASTA/summary consistency
- Scaffold exclusion
- UTR removal (CDS extraction)

Usage:
    python validate_fasta_output.py output.fasta summary.tsv [--min_length 26]
"""

import argparse
import sys
import logging
from collections import Counter

import pandas as pd

# Setup basic logger structure, config will be updated in main
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s | %(message)s',
    datefmt='%H:%M:%S'
)

logger = logging.getLogger(__name__)

def parse_fasta(fasta_path):
    """Parse FASTA file and return list of records."""
    records = []
    with open(fasta_path, 'r') as f:
        header = None
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header and seq:
                    records.append({'header': header, 'seq': seq})
                header = line
                seq = ''
            else:
                seq += line
        if header and seq:
            records.append({'header': header, 'seq': seq})
    return records


def validate_fasta(fasta_path, summary_path, min_length=26):
    """Run all validation checks and return success status."""
    
    fasta_records = parse_fasta(fasta_path)
    summary_df = pd.read_csv(summary_path, sep='\t', low_memory=False)
    all_passed = True
    
    logger.info(f"=== FASTA Validation ({len(fasta_records):,} records) ===")

    # 1. Header Format
    bad_headers, biotypes, transcript_ids = [], [], []
    for r in fasta_records:
        parts = r['header'][1:].split('|')
        if len(parts) != 5:
            bad_headers.append(r['header'])
        else:
            biotypes.append(parts[4])
            transcript_ids.append(parts[2])
    
    if bad_headers:
        logger.error(f"[FAIL] Header format: {len(bad_headers)} invalid (e.g. {bad_headers[0]})")
        all_passed = False
    else:
        logger.info("[PASS] Header format")

    # 2. Biotype
    biotype_counts = Counter(biotypes)
    invalid_biotypes = [b for b in biotypes if b not in ('protein_coding', 'lncRNA')]
    if invalid_biotypes:
        logger.error(f"[FAIL] Biotypes: {len(invalid_biotypes)} invalid")
        all_passed = False
    else:
        logger.info(f"[PASS] Biotypes: {dict(biotype_counts)}")

    # 3. Sequence validity
    invalid_seqs = [(r['header'], set(r['seq']) - set('ACTGN')) 
                    for r in fasta_records if not all(c in 'ACTGN' for c in r['seq'])]
    if invalid_seqs:
        logger.error(f"[FAIL] Sequence chars: {len(invalid_seqs)} invalid")
        all_passed = False
    else:
        logger.info("[PASS] Sequence chars (ACTGN only)")

    # 4. Sequence length
    lengths = [len(r['seq']) for r in fasta_records]
    short_seqs = [l for l in lengths if l < min_length]
    if short_seqs:
        logger.error(f"[FAIL] Length: {len(short_seqs)} < {min_length}bp")
        all_passed = False
    else:
        logger.info(f"[PASS] Length: min={min(lengths)}, max={max(lengths):,}, med={sorted(lengths)[len(lengths)//2]:,}")

    # 5. Duplicates
    dup_ids = [tid for tid, count in Counter(transcript_ids).items() if count > 1]
    if dup_ids:
        logger.error(f"[FAIL] Duplicates: {len(dup_ids)} (e.g. {dup_ids[0]})")
        all_passed = False
    else:
        logger.info("[PASS] No duplicates")

    # 6. FASTA-Summary consistency
    included = summary_df[summary_df['included_in_output_fasta'] == True]
    fasta_tx_set, summary_tx_set = set(transcript_ids), set(included['transcript_id'])
    diff1, diff2 = fasta_tx_set - summary_tx_set, summary_tx_set - fasta_tx_set
    if diff1 or diff2:
        logger.error(f"[FAIL] FASTA/summary mismatch: +{len(diff1)}/-{len(diff2)}")
        all_passed = False
    else:
        logger.info(f"[PASS] FASTA/summary match ({len(fasta_tx_set):,} transcripts)")

    # 7. Scaffold exclusion
    scaffolds = set(summary_df[summary_df['exclusion_reason'] == 'scaffold_patch_haplotype']['transcript_id'])
    scaffolds_in_fasta = fasta_tx_set & scaffolds
    if scaffolds_in_fasta:
        logger.error(f"[FAIL] Scaffolds: {len(scaffolds_in_fasta)} in output")
        all_passed = False
    else:
        logger.info(f"[PASS] Scaffolds excluded ({len(scaffolds):,})")

    # 8. CDS extraction
    included_pc = summary_df[(summary_df['included_in_output_fasta']) & (summary_df['biotype'] == 'protein_coding')]
    if len(included_pc) > 0 and 'cds_length' in summary_df.columns:
        pc_in_fasta = {r['header'].split('|')[2]: len(r['seq']) 
                       for r in fasta_records if r['header'].split('|')[4] == 'protein_coding'}
        mismatches = [(tid, pc_in_fasta[tid], int(row['cds_length'])) 
                      for _, row in included_pc.head(100).iterrows() 
                      if row['transcript_id'] in pc_in_fasta and pc_in_fasta[row['transcript_id']] != int(row['cds_length'])]
        if mismatches:
            logger.error(f"[FAIL] CDS length: {len(mismatches)} mismatches")
            all_passed = False
        else:
            utrs_removed = (included_pc['cds_length'] < included_pc['transcript_length']).sum()
            logger.info(f"[PASS] CDS extraction ({utrs_removed:,} UTRs removed)")
    else:
        logger.info("[SKIP] CDS check (no data)")

    # 9. Representative Transcript Checks
    if 'is_representative' in summary_df.columns:
        reps = summary_df[summary_df['is_representative'] == True]
        logger.info(f"--- Representative Validation ({len(reps):,}) ---")
        
        # Check 9a: Must have Ensembl_canonical tag
        if 'tag' in summary_df.columns:
            non_canonical_reps = reps[~reps['tag'].astype(str).str.contains('Ensembl_canonical', na=False)]
            if not non_canonical_reps.empty:
                logger.error(f"[FAIL] Representatives without Ensembl_canonical tag: {len(non_canonical_reps)}")
                logger.error(f"       Examples: {non_canonical_reps['transcript_id'].head().tolist()}")
                all_passed = False
            else:
                logger.info("[PASS] All representatives have Ensembl_canonical tag")
        
        # Check 9b: Must have full length preserved (no CDS trimming)
        # We expect cds_length == transcript_length for representatives
        if 'cds_length' in summary_df.columns and 'transcript_length' in summary_df.columns:
            trimmed_reps = reps[reps['cds_length'] != reps['transcript_length']]
            if not trimmed_reps.empty:
                logger.error(f"[FAIL] Representatives with trimmed sequences: {len(trimmed_reps)}")
                logger.error(f"       Examples: {trimmed_reps['transcript_id'].head().tolist()}")
                all_passed = False
            else:
                logger.info("[PASS] All Ensembl_canonical representatives preserved full length (UTRs included)")

    # Summary stats
    logger.info(f"--- Stats ---")
    if 'exclusion_reason' in summary_df.columns:
        excl = summary_df[summary_df['exclusion_reason'].notna()]['exclusion_reason'].value_counts()
        logger.info(f"Excluded: {', '.join(f'{r}={c:,}' for r, c in excl.items())}")
    logger.info(f"Output: {', '.join(f'{b}={c:,}' for b, c in biotype_counts.items())}")

    logger.info(f"{'PASSED' if all_passed else 'FAILED'}")
    return all_passed


def main():
    parser = argparse.ArgumentParser(
        description='Validate FASTA reformatter output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Examples:
                python validate_fasta_output.py output.fasta summary.tsv
                python validate_fasta_output.py output.fasta summary.tsv --min_length 100
        """
    )
    
    parser.add_argument('fasta', help='Output FASTA file to validate')
    parser.add_argument('summary', help='Summary TSV file')
    parser.add_argument('--min_length', '-m', type=int, default=26,
                        help='Minimum sequence length (default: 26)')
    parser.add_argument('--log', '-l', default='validation.log',
                        help='Log file path (default: validation.log)')
    
    args = parser.parse_args()
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(args.log),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    success = validate_fasta(args.fasta, args.summary, args.min_length)
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()