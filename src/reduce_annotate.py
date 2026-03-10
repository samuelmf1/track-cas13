#!/usr/bin/env python3
import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv', required=True)
    parser.add_argument('--exons', required=True)
    parser.add_argument('--cds', required=True)
    parser.add_argument('--summary', required=True)
    parser.add_argument('--expr', required=True)
    parser.add_argument('--out_prefix', required=True)
    args = parser.parse_args()

    # 1. LOAD DATA (Read once, use many)
    df = pd.read_csv(args.csv, low_memory=False)
    
    # Load lookups and set indexes for O(1) access
    summary_df = pd.read_csv(args.summary, sep='\t').set_index('transcript_id')
    expr_df = pd.read_csv(args.expr, sep='\t').set_index('transcript_id')
    cds_df = pd.read_csv(args.cds, sep='\t')
    cds_dict = cds_df[cds_df['region'] == 'CDS'].set_index('transcript_id')[['start', 'end']].to_dict('index')
    
    # Pre-process Exons into a dictionary of lists for fast iteration per transcript
    exon_raw = pd.read_csv(args.exons, sep='\t')
    exon_lookup = {}
    for tid, group in exon_raw.groupby('transcript_id'):
        exon_lookup[tid] = list(zip(group['exon_id'], group['exon_sequence'].str.upper()))

    # 2. FAST EXON MATCHING
    def find_exon(row):
        tid = row['Symbol']
        seq = row['Target Sequence'].upper()
        if tid in exon_lookup:
            for eid, eseq in exon_lookup[tid]:
                if seq in eseq: return eid
        return None

    df['exon_id'] = df.apply(find_exon, axis=1)
    
    # Create the "a1" output (Exon annotated)
    # df.to_csv(f"{args.out_prefix}.a1.csv", index=False)

    # 3. FILTER & REGION ANNOTATION (CDS/UTR)
    # Filter "inadmissible" (equivalent to your grep -v)
    df = df[~df['exon_id'].str.contains('inadmissible', case=False, na=False)].copy()

    def get_region(row, min_cds_overlap=23):
        """
        Classifies the guide target region. 
        min_cds_overlap defaults to 23 for strict 100% CDS classification.
        Change to 12 for majority (>50%) classification.
        """
        tid = row['Symbol']
        if tid not in summary_df.index: 
            return "Unknown"
        
        info = summary_df.loc[tid]
        idx = row['Contig_Idx']
        
        # Adjust index for representative transcripts 
        if info['biotype'] == 'lncRNA': 
            return "lncRNA"
        # if info['biotype'] == 'protein_coding' and not info['is_representative']:
        #     return "CDS"
        if info['biotype'] == 'protein_coding':# and info['is_representative']:
            idx += 3
        
        if tid not in cds_dict: 
            return info['biotype']
        
        c = cds_dict[tid]
        
        # Guide boundaries (assuming a 23nt spacer)
        g_start, g_end = idx, idx + 23
        
        # Calculate exact nucleotide overlap with the CDS
        overlap = max(0, min(g_end, c['end']) - max(g_start, c['start']))
        
        # Check against our defined stringency threshold
        if overlap >= min_cds_overlap: 
            return "CDS"
        
        # If it fails the CDS threshold, bin it into the appropriate UTR.
        # If the guide starts before the CDS, it's anchored in the 5'UTR.
        # Otherwise, it's anchored in (or entirely within) the 3'UTR.
        return "5'UTR" if g_start < c['start'] else "3'UTR"
    
    df['region'] = df.apply(get_region, axis=1)
    # df.to_csv(f"{args.out_prefix}.a2.csv", index=False)

    # 4. TAGS & EXPRESSION (Vectorized Mapping)
    df['tags'] = df['Symbol'].map(
        summary_df['tag']
        .str.replace(';', ',', regex=False)  # Normalize delimiters to commas
        .str.replace(',', '|', regex=False)  # Convert commas to pipes
        .str.strip()
        .to_dict()
    )
    df['median_expr'] = df['Symbol'].map(expr_df['median_expr'].to_dict()).fillna(-1)
    df['median_expr_norm'] = df['Symbol'].map(expr_df['median_expr_norm'].to_dict()).fillna(-1)
    df['is_expressed'] = df['Symbol'].map(expr_df['is_expressed'].to_dict()).fillna(False)
    df['pce'] = df['Symbol'].map(expr_df['pce'].to_dict()).fillna(-1)
    df['n_expressed'] = df['Symbol'].map(expr_df['n_expressed'].to_dict()).fillna(-1)

    df.to_csv(f"{args.out_prefix}.annots.csv", index=False)

if __name__ == "__main__":
    main()