#!/usr/bin/env python3
import sys
import pandas as pd


# --- GTF Parsing Functions ---

LOW_CONF_TAGS = ['low_conf', 'partial', 'fragment', 'artifact', 'cds_start_NF', 'cds_end_NF']
HIGH_CONF_TAGS = ['basic', 'CCDS']

def parse_gtf_attributes(attr_string: str) -> pd.Series:
    """Parses GTF attribute string into key-value pairs."""
    data = {
        'gene_id': None,
        'transcript_id': None,
        'gene_name': None,
        'gene_type': None,
        'transcript_type': None,
        'transcript_support_level': None,
        'tag': [],
    }

    for item in attr_string.split('; '):
        if not item:
            continue
        try:
            key, value = item.split(' ', 1)
            value = value.strip().replace('"', '').rstrip(';')
            if key == 'tag':
                data['tag'].append(value)
            elif key in data:
                data[key] = value
        except ValueError:
            continue
    
    data['tag'] = ','.join(data['tag']) if data['tag'] else ''
    return pd.Series(data)

def has_no_low_conf_tag(tags_string: str) -> bool:
    """Returns True if NONE of the low-confidence tags are present."""
    if not tags_string: 
        return True
    tag_list = [t.strip() for t in tags_string.split(',')]
    return not any(tag in tag_list for tag in LOW_CONF_TAGS)


def has_basic_or_ccds_tag(tags_string: str) -> bool:
    """Returns True if EITHER 'basic' OR 'CCDS' tag is present."""
    if not tags_string: 
        return False
    tag_list = [t.strip() for t in tags_string.split(',')]
    return any(tag in tag_list for tag in HIGH_CONF_TAGS)


def load_gtf_mapping(gtf_path: str) -> pd.DataFrame:
    """
    Load and parse GTF file to create transcript-to-gene mapping with annotations.
    """
    print(f"Loading GTF file: {gtf_path}...", file=sys.stderr)
    gtf_cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    
    gtf_df = pd.read_csv(
        gtf_path, 
        sep='\t', 
        comment='#', 
        header=None, 
        names=gtf_cols, 
        low_memory=False,
        compression='infer',
        usecols=['feature', 'attribute']
    )

    df_transcripts = gtf_df[gtf_df['feature'] == 'transcript'].copy()
    print(f"  Transcript features in GTF: {len(df_transcripts):,}", file=sys.stderr)
    
    attr_df = df_transcripts['attribute'].apply(parse_gtf_attributes)
    df_transcripts = pd.concat([df_transcripts.drop('attribute', axis=1), attr_df], axis=1)
    df_transcripts['biotype'] = df_transcripts['transcript_type']
    
    mapping_df = df_transcripts[[
        'transcript_id', 
        'gene_id', 
        'gene_name',
        'biotype', 
        'transcript_support_level',
        'tag'
    ]].drop_duplicates(subset=['transcript_id']).reset_index(drop=True)

    mapping_df['tsl_numeric'] = pd.to_numeric(mapping_df['transcript_support_level'], errors='coerce')
    
    # Identify representative transcript per gene (solely based on Ensembl_canonical tag)
    # We want to prioritize transcripts that have the 'Ensembl_canonical' tag.
    # If multiple have it (unlikely), valid tie-breakers are internal consistency (id).
    
    def is_ensembl_canonical(tag_str):
        if pd.isna(tag_str) or not isinstance(tag_str, str):
            return False
        return 'Ensembl_canonical' in [t.strip() for t in tag_str.split(',')]

    mapping_df['is_ensembl_canonical'] = mapping_df['tag'].apply(is_ensembl_canonical)
    
    # Sort to bring Ensembl_canonical to the top
    mapping_df.sort_values(
        by=['gene_id', 'is_ensembl_canonical', 'transcript_id'], 
        ascending=[True, False, True], 
        inplace=True
    )
    
    # Mark the first one for each gene as representative IF it is canonical
    mapping_df['is_representative'] = False
    
    # We only want to set is_representative=True if the top one IS actually canonical
    # The sort puts True before False. So if a gene has a canonical tx, it's at the top.
    # If it doesn't, a non-canonical is at the top.
    
    # Get the top entry for each gene
    top_entries = mapping_df.groupby('gene_id').head(1)
    
    # Filter for those that are actually canonical
    canonical_representatives = top_entries[top_entries['is_ensembl_canonical']].index
    
    mapping_df.loc[canonical_representatives, 'is_representative'] = True

    # Independent filtering flags
    # tsl_ok: TSL 1-4 (good transcript support level)
    mapping_df['tsl_ok'] = (mapping_df['tsl_numeric'] <= 4)
    # no_low_conf_tags: no low-confidence tags present
    mapping_df['no_low_conf_tags'] = mapping_df['tag'].apply(has_no_low_conf_tag)
    # is_basic_or_ccds: has basic or CCDS tag
    mapping_df['is_basic_or_ccds'] = mapping_df['tag'].apply(has_basic_or_ccds_tag)
    
    mapping_df = mapping_df[mapping_df['biotype'].notna()]
    
    print(f"  Parsed transcripts: {len(mapping_df):,}", file=sys.stderr)
    print(f"    Protein-coding: {(mapping_df['biotype'] == 'protein_coding').sum():,}", file=sys.stderr)
    print(f"    lncRNA: {(mapping_df['biotype'] == 'lncRNA').sum():,}", file=sys.stderr)
    print(f"    Representatives (Ensembl_canonical): {mapping_df['is_representative'].sum():,}", file=sys.stderr)

    return mapping_df.drop(columns=['tsl_numeric', 'is_ensembl_canonical']).rename(columns={'biotype': 'gene_type'})


# --- FASTA Processing Functions ---

def parse_header(header):
    """Extracts ENSG|ENST|GENE|BIOTYPE from GENCODE header"""
    parts = header[1:].split('|')
    enst = parts[0] if len(parts) > 0 else 'unknown'
    ensg = parts[1] if len(parts) > 1 else 'unknown'
    gene = parts[5] if len(parts) > 5 else 'unknown'
    biotype = parts[7] if len(parts) > 7 else (parts[6] if len(parts) > 6 else 'unknown')
    return f">{ensg}|{enst}|{gene}|{biotype}"


def is_valid_dna(seq):
    """Check if sequence contains only A, C, T, G, N (case-insensitive)"""
    return all(c.upper() in 'ACTGN' for c in seq)


def write_fasta_record(outfile, header, sequence):
    """Write a FASTA record with exactly one header line and one sequence line."""
    outfile.write(header.rstrip() + '\n')
    outfile.write(sequence.replace('\n', '').replace('\r', '') + '\n')


def process_transcript(header, seq_lines, transcript_counter, outfile, min_length, 
                       abundant_set, biotype_set, summary_data, cds_dict, 
                       gtf_transcript_set, gtf_basic_ccds_set, representative_set, longest_utr_set):
    if not header or not seq_lines:
        return transcript_counter

    full_seq = ''.join(seq_lines)
    seq_length = len(full_seq)
    
    # Parse original header fields
    try:
        parts = header[1:].split('|')
        ensg = parts[1] if len(parts) > 1 else 'unknown'
        enst = parts[0] if len(parts) > 0 else 'unknown'
        gene = parts[5] if len(parts) > 5 else 'unknown'
        biotype = parts[7] if len(parts) > 7 else (parts[6] if len(parts) > 6 else 'unknown')
    except Exception as e:
        print(f"ERROR: Failed to parse header {header}: {e}", file=sys.stderr)
        summary_data.append({
            'gene_id': 'parse_error',
            'transcript_id': 'parse_error',
            'biotype': 'unknown',
            'transcript_length': seq_length,
            'cds_length': 0,
            'included_in_output_fasta': False,
            'exclusion_reason': 'header_parse_error',
            'is_representative': False,
            'is_longest_utr': False
        })
        return transcript_counter

    # Check if this transcript is the designated representative or longest UTR for its gene
    is_rep = representative_set is not None and enst in representative_set
    is_longest_utr = longest_utr_set is not None and enst in longest_utr_set
    
    # Preservation flag: keep UTRs for representative OR longest UTR transcripts
    keep_utrs = is_rep or is_longest_utr

    # Initialize exclusion tracking
    included = True
    exclusion_reason = None
    
    # Check biotype first (Representatives MUST still match biotype)
    if biotype not in biotype_set:
        included = False
        exclusion_reason = 'wrong_biotype'
    # Check if transcript is in GTF (not scaffold/patch/haplotype) - only for pc/lncRNA
    # Representatives MUST be in GTF (since we picked them from GTF)
    elif gtf_transcript_set is not None and enst not in gtf_transcript_set:
        included = False
        exclusion_reason = 'scaffold_patch_haplotype'
    # Check if transcript has basic or CCDS tag
    # BYPASS: Preservation transcripts bypass this filter
    elif not keep_utrs and biotype == 'protein_coding' and gtf_basic_ccds_set is not None and enst not in gtf_basic_ccds_set:
        included = False
        exclusion_reason = 'not_basic_ccds_pc_tx'
    # Check sequence length
    elif seq_length < min_length:
        # print(f"WARNING: Skipping short transcript {header} (len={seq_length})", file=sys.stderr)
        included = False
        exclusion_reason = f'too_short (tx len<{min_length})'
    # Check for non-DNA characters
    elif not is_valid_dna(full_seq):
        print(f"WARNING: Skipping transcript with non-DNA chars: {header}", file=sys.stderr)
        included = False
        exclusion_reason = 'non_DNA_characters'
    # Check abundance
    # BYPASS: Preservation transcripts bypass abundance filter
    elif not keep_utrs and abundant_set is not None and enst not in abundant_set:
        included = False
        exclusion_reason = 'not_abundant'
    
    # Keep full sequence (includes UTRs for all transcripts at this stage)
    cds_len = seq_length
    cds_seq = full_seq

    # Add to summary data
    summary_data.append({
        'gene_id': ensg,
        'transcript_id': enst,
        'biotype': biotype,
        'transcript_length': seq_length,
        'cds_length': len(cds_seq), # Actual written length
        'included_in_output_fasta': included,
        'exclusion_reason': exclusion_reason if not included else None,
        'is_representative': is_rep,
        'is_longest_utr': is_longest_utr
    })
    
    # Write to FASTA if included
    if included:
        transcript_counter += 1
        new_header = f">{transcript_counter}|{ensg}|{enst}|{gene}|{biotype}"
        write_fasta_record(outfile, new_header, cds_seq)

    return transcript_counter


def main(input_path, output_path, cds_boundaries_path=None, abundance_path=None, 
        abundance_column='is_abundant', min_length=30, summary_path=None, gtf_path=None):
    
    # Load GTF annotations if provided
    gtf_df = None
    gtf_transcript_set = None
    gtf_basic_ccds_set = None
    representative_set = None
    longest_utr_set = None
    
    if gtf_path:
        gtf_df = load_gtf_mapping(gtf_path)
        gtf_transcript_set = set(gtf_df['transcript_id'])
        # Create set of transcripts with basic or CCDS tag
        gtf_basic_ccds_set = set(gtf_df[gtf_df['is_basic_or_ccds']]['transcript_id'])
        # Create set of representative transcripts for UTR preservation
        representative_set = set(gtf_df[gtf_df['is_representative']]['transcript_id'])
        
        print(f"  Will filter out transcripts not in GTF (scaffold/patch/haplotype)", file=sys.stderr)
        print(f"  Will filter out transcripts without basic/CCDS tag ({len(gtf_basic_ccds_set):,} have tag)", file=sys.stderr)
        print(f"  Identified {len(representative_set):,} representative transcripts (one per gene) for UTR preservation.", file=sys.stderr)
    else:
        print("No GTF file provided - scaffold/patch/haplotype and basic/CCDS filtering disabled.", file=sys.stderr)

    # Load CDS boundaries and identify longest UTR transcripts
    cds_dict = {}
    if cds_boundaries_path:
        print(f"Loading CDS boundaries from {cds_boundaries_path}...", file=sys.stderr)
        cds_boundaries_df = pd.read_table(cds_boundaries_path)
        
        # Calculate UTR lengths for each transcript
        utr_data = cds_boundaries_df[cds_boundaries_df['region'].isin(['5UTR', '3UTR'])].copy()
        utr_data['length'] = utr_data['end'] - utr_data['start'] + 1
        utr_lengths = utr_data.groupby('transcript_id')['length'].sum().reset_index()
        
        # Merge with GTF to get gene association
        if gtf_df is not None:
            tx_to_gene = gtf_df[['transcript_id', 'gene_id']].drop_duplicates()
            utr_lengths = utr_lengths.merge(tx_to_gene, on='transcript_id', how='left')
            
            # Identify transcript with max UTR length per gene
            longest_utr_idx = utr_lengths.sort_values(['gene_id', 'length', 'transcript_id'], 
                                                     ascending=[True, False, True]).groupby('gene_id').head(1).index
            longest_utr_set = set(utr_lengths.loc[longest_utr_idx, 'transcript_id'])
            print(f"  Identified {len(longest_utr_set):,} transcripts with the longest UTR for their respective genes.", file=sys.stderr)
        else:
            print("WARNING: Cannot identify longest UTR per gene without GTF; longest UTR preservation disabled.", file=sys.stderr)

        # Create dictionary for CDS trimming
        cds_dict = cds_boundaries_df[cds_boundaries_df['region'] == 'CDS'].copy().set_index('transcript_id').to_dict(orient='index')
    else:
        print("No CDS boundaries file provided - keeping full length of transcripts that meet other criteria.", file=sys.stderr)

    # Load abundance data if provided
    abundant_set = None
    skip_abundance_values = {'NA', 'na', 'None', 'none', 'skip', 'SKIP', ''}
    
    if abundance_path is not None:
        if abundance_column in skip_abundance_values:
            print(f"Abundance file provided but column set to '{abundance_column}' - skipping abundance filtering.", file=sys.stderr)
        else:
            try:
                print(f"Loading abundance data from {abundance_path}...", file=sys.stderr)
                abundance_df = pd.read_table(abundance_path)
                
                if abundance_column not in abundance_df.columns:
                    available_columns = ', '.join(abundance_df.columns)
                    raise ValueError(f"Abundance TSV must contain '{abundance_column}' column. Available columns: {available_columns}")
                
                if 'transcript_id' not in abundance_df.columns:
                    raise ValueError("Abundance TSV must contain 'transcript_id' column")
                
                abundant_df = abundance_df[abundance_df[abundance_column]]
                abundant_set = set(abundant_df['transcript_id'])
                print(f"Found {len(abundant_set)} abundant transcripts to filter by.", file=sys.stderr)
            except Exception as e:
                print(f"ERROR: Failed to load abundance file {abundance_path}: {e}", file=sys.stderr)
                sys.exit(1)
    else:
        print("No abundance file provided - keeping all transcripts that meet other criteria.", file=sys.stderr)

    biotype_set = {'lncRNA', 'protein_coding'}

    transcript_counter = 0
    current_header = None
    current_seq_lines = []
    total_transcripts_seen = 0
    summary_data = []

    try:
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for raw_line in infile:
                line = raw_line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    if current_header:
                        total_transcripts_seen += 1
                        transcript_counter = process_transcript(
                            current_header, current_seq_lines,
                            transcript_counter, outfile, min_length,
                            abundant_set, biotype_set, summary_data, cds_dict, 
                            gtf_transcript_set, gtf_basic_ccds_set, representative_set, longest_utr_set
                        )
                    current_header = line
                    current_seq_lines = []
                else:
                    current_seq_lines.append(line.upper())

            # Process last transcript
            if current_header:
                total_transcripts_seen += 1
                transcript_counter = process_transcript(
                    current_header, current_seq_lines,
                    transcript_counter, outfile, min_length,
                    abundant_set, biotype_set, summary_data, cds_dict, 
                    gtf_transcript_set, gtf_basic_ccds_set, representative_set, longest_utr_set
                )

        print(f"Success: Processed {transcript_counter}/{total_transcripts_seen} transcripts -> {output_path}", file=sys.stderr)
        
        # Create summary DataFrame and merge with GTF annotations
        if summary_path:
            summary_df = pd.DataFrame(summary_data)
            
            # Merge with GTF annotations if available
            if gtf_df is not None:
                # Rename columns to avoid conflicts
                # Rename columns to avoid conflicts
                gtf_merge_cols = ['transcript_id', 'gene_name', 'gene_type', 'transcript_support_level', 
                                  'tag', 'tsl_ok', 'no_low_conf_tags', 'is_basic_or_ccds', 'is_representative']
                
                # We'll calculate is_longest_utr status separately from the longest_utr_set
                
                gtf_for_merge = gtf_df[gtf_merge_cols].copy()
                
                # Merge - left join to keep all FASTA transcripts
                # Drop overlapping flag columns from summary_df to overwrite with GTF or re-calculated versions
                cols_to_drop = [c for c in ['is_representative', 'is_longest_utr'] if c in summary_df.columns]
                summary_df = summary_df.drop(columns=cols_to_drop)
                     
                summary_df = summary_df.merge(gtf_for_merge, on='transcript_id', how='left')
                
                # Re-add is_longest_utr flag
                summary_df['is_longest_utr'] = False
                if longest_utr_set is not None:
                     summary_df.loc[summary_df['transcript_id'].isin(longest_utr_set), 'is_longest_utr'] = True
                
                # Fill NaN for transcripts not in GTF (scaffolds)
                summary_df['gene_type'] = summary_df['gene_type'].fillna(summary_df['biotype'])
                summary_df['is_representative'] = summary_df['is_representative'].fillna(False)
            
            summary_df.to_csv(summary_path, sep='\t', index=False)
            print(f"Summary written to {summary_path}", file=sys.stderr)
            
            # Print summary statistics
            total_parsed = len(summary_df)
            total_included = summary_df['included_in_output_fasta'].sum()
            print(f"\nSummary Statistics:", file=sys.stderr)
            print(f"  Total transcripts parsed: {total_parsed}", file=sys.stderr)
            print(f"  Transcripts included in output: {total_included}", file=sys.stderr)
            print(f"  Transcripts excluded: {total_parsed - total_included}", file=sys.stderr)
            
            if 'exclusion_reason' in summary_df.columns:
                exclusion_counts = summary_df[summary_df['exclusion_reason'].notna()]['exclusion_reason'].value_counts()
                if not exclusion_counts.empty:
                    print(f"\nExclusion reasons:", file=sys.stderr)
                    for reason, count in exclusion_counts.items():
                        print(f"    {reason}: {count}", file=sys.stderr)
        
        # Report filtering summary
        filters_applied = ["biotype (lncRNA/protein_coding)", f"min_length ({min_length})"]
        if gtf_transcript_set is not None:
            filters_applied.append("GTF primary assembly filter")
        if gtf_basic_ccds_set is not None:
            filters_applied.append("basic/CCDS tag filter")
        if abundant_set is not None:
            filters_applied.append("abundance filter")
        print(f"\nFiltering applied: {', '.join(filters_applied)}", file=sys.stderr)

    except Exception as e:
        print(f"ERROR: Failed during processing: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Transform FASTA file with GTF annotation and optional abundance filtering',
        usage='%(prog)s <input.fa> <output.fasta> [options]'
    )

    parser.add_argument('input_fa', help='Input FASTA file')
    parser.add_argument('output_fasta', help='Output FASTA file')
    
    parser.add_argument('--gtf', '-g', default=None,
                        help='GTF annotation file to filter scaffold/patch/haplotype transcripts')
    
    parser.add_argument('--cds_boundaries', '-b', default=None,
                        help='TSV file containing CDS boundaries (optional)')

    parser.add_argument('--abundance_tsv', '-a', default=None,
                        help='TSV file containing abundance data (optional)')
    parser.add_argument('--abundance_column', '-c', default='is_abundant',
                        help='Name of the abundance column to use (default: is_abundant). Use "NA" to skip filtering')
    
    parser.add_argument('--min_length', '-m', type=int, default=26,
                        help='Minimum sequence length (default: 26)')
    
    parser.add_argument('--summary', '-s', default=None,
                        help='Output TSV file for summary of all parsed transcripts (optional)')

    args = parser.parse_args()

    if args.input_fa == args.output_fasta:
        print("ERROR: Input and output paths must be different.", file=sys.stderr)
        sys.exit(1)

    main(args.input_fa, args.output_fasta, args.cds_boundaries, args.abundance_tsv, 
         args.abundance_column, args.min_length, args.summary, args.gtf)


"""
GENCODE FASTA reformatter with GTF annotation integration:
- One sequence per line
- Skips lines with non-DNA characters
- Skips sequences shorter than minimum length
- Keeps only lncRNA and protein_coding transcripts
- Filters out scaffold/patch/haplotype transcripts (not in primary GTF)
- Filters out transcripts without basic or CCDS tag
- OPTIONALLY keeps only abundant transcripts (from abundance TSV if provided)
  - Can skip abundance filtering by setting column to "NA" even with file provided
- Clean header: >1|ENSG|ENST|GENE|BIOTYPE
- Outputs summary TSV with all parsed transcripts, inclusion status, and GTF annotations
"""