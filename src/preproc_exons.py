#!/usr/bin/env python3
"""
Script to extract all exons from a GTF file and their sequences from a FASTA file.
Also extracts junction sequences with seven types per junction:
  1. Admissible: -19:+19 (balanced, last 19bp of exon N + first 19bp of exon N+1)
  2-7. Inadmissible (asymmetric, guides barely cross the junction boundary):
       -22:1, -21:2, -20:3, -3:20, -2:21, -1:22

       ** we have changed this to be more stringent since 
       ** https://academic.oup.com/nar/article/54/1/gkaf1447/8417327?login=false

Plus special handling for short middle exons (< 19bp):
  - Creates 3-exon sequence: exon1[-19:] + exon2[full] + exon3[:+19]
  - Plus inadmissible 2-exon junctions with the short exon

Inadmissible junctions are marked with "_inadmissible_LEFT:RIGHT" suffix and is_admissible=False.

Outputs: exon_id, transcript_id, gene_id, exon_number, strand, sequence, length, is_admissible
"""

import re
import os
import sys
import argparse
import time
from collections import defaultdict
from typing import Dict, List, Tuple, Set


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract exon and junction sequences from GTF and FASTA files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -g annotations.gtf -f transcripts.fa -o exon_sequences.tsv
  %(prog)s --gtf annotations.gtf --fasta transcripts.fa --output results.tsv --verbose
        """
    )
    
    # Required arguments
    parser.add_argument(
        "-g", "--gtf",
        required=True,
        help="Input GTF file with exon annotations"
    )
    parser.add_argument(
        "-f", "--fasta",
        required=True,
        help="Input FASTA file with transcript sequences"
    )
    
    # Optional arguments
    parser.add_argument(
        "-o", "--output",
        default="exon_sequences.tsv",
        help="Output TSV file (default: exon_sequences.tsv)"
    )
    parser.add_argument(
        "--exons-only",
        action="store_true",
        help="Extract only exons, skip junction sequences"
    )
    parser.add_argument(
        "--junctions-only",
        action="store_true",
        help="Extract only junctions, skip individual exon sequences"
    )
    parser.add_argument(
        "--skip-inadmissible",
        action="store_true",
        help="Skip inadmissible (asymmetric) junction sequences"
    )
    parser.add_argument(
        "--short-exon-threshold",
        type=int,
        default=19,
        help="Threshold for short middle exon handling (default: 19)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print verbose progress information"
    )
    parser.add_argument(
        "--max-warnings",
        type=int,
        default=10,
        help="Maximum number of warnings to print (default: 10)"
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.exons_only and args.junctions_only:
        parser.error("Cannot use both --exons-only and --junctions-only")
    
    # Check input files exist
    if not os.path.isfile(args.gtf):
        parser.error(f"GTF file not found: {args.gtf}")
    if not os.path.isfile(args.fasta):
        parser.error(f"FASTA file not found: {args.fasta}")
    
    # Check output directory is writable
    output_dir = os.path.dirname(args.output) or "."
    if not os.access(output_dir, os.W_OK):
        parser.error(f"Output directory not writable: {output_dir}")
    
    return args


def parse_gtf_attribute(attribute_string: str, key: str) -> str:
    """Extract a specific attribute value from GTF attribute string."""
    pattern = f'{key} "([^"]+)"'
    match = re.search(pattern, attribute_string)
    return match.group(1) if match else None


def load_fasta_sequences(fasta_file: str, verbose: bool = False) -> Dict[str, str]:
    """Load transcript sequences from FASTA file into dictionary.
    
    Supports two header formats:
    1. Original GENCODE: >ENST00000832824.1|ENSG00000290825.2|...|lncRNA|
       - transcript_id is first field (position 0)
    2. Reformatted subset: >1|ENSG00000238009.6|ENST00000466430.5|RP11-34P13.7|lncRNA
       - transcript_id is third field (position 2)
    """
    sequences = {}
    current_transcript = None
    current_seq = []
    
    if verbose:
        print("  Reading FASTA file...")
    
    # Auto-detect format from first header
    header_format = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Save previous sequence
                if current_transcript:
                    sequences[current_transcript] = ''.join(current_seq)
                
                header = line[1:].rstrip()
                parts = header.split('|')
                
                # Auto-detect format on first header
                if header_format is None:
                    if len(parts) >= 1 and parts[0].startswith('ENST'):
                        header_format = 'gencode_original'
                        if verbose:
                            print(f"  Detected GENCODE original header format")
                    else:
                        header_format = 'reformatted'
                        if verbose:
                            print(f"  Detected reformatted header format")
                
                # Parse based on detected format
                if header_format == 'gencode_original':
                    # Original GENCODE: >ENST00000832824.1|ENSG00000290825.2|...
                    current_transcript = parts[0] if len(parts) >= 1 else header.split()[0]
                else:
                    # Reformatted: >1|ENSG...|ENST...|...
                    current_transcript = parts[2] if len(parts) >= 3 else header.split()[0]
                
                current_seq = []
            else:
                current_seq.append(line.rstrip())
        
        # Save last sequence
        if current_transcript:
            sequences[current_transcript] = ''.join(current_seq)
    
    return sequences


def parse_gtf_exons(gtf_file: str, verbose: bool = False) -> Tuple[Dict[str, List[Tuple[int, int, str, int]]], Dict[str, str], Dict[str, str]]:
    """
    Parse GTF file and extract exons for each transcript.
    Returns: transcript_exons (with exon_id and exon_number), transcript_strand, transcript_gene
    """
    transcript_exons = defaultdict(list)
    transcript_strand = {}
    transcript_gene = {}
    
    # Pre-compile regex patterns
    transcript_pattern = re.compile(r'transcript_id "([^"]+)"')
    gene_pattern = re.compile(r'gene_id "([^"]+)"')
    exon_id_pattern = re.compile(r'exon_id "([^"]+)"')
    exon_number_pattern = re.compile(r'exon_number (\d+)')
    
    if verbose:
        print("  Reading GTF file...")
    
    line_count = 0
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            line_count += 1
            if verbose and line_count % 100000 == 0:
                print(f"    Processed {line_count:,} lines...", end='\r')
            
            fields = line.rstrip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            
            # Only process exons
            if feature_type != 'exon':
                continue
            
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            # Extract transcript_id
            match = transcript_pattern.search(attributes)
            if not match:
                continue
            transcript_id = match.group(1)
            
            # Extract gene_id
            gene_match = gene_pattern.search(attributes)
            gene_id = gene_match.group(1) if gene_match else "NA"
            
            # Extract exon_id (ENSE number)
            exon_id_match = exon_id_pattern.search(attributes)
            exon_id = exon_id_match.group(1) if exon_id_match else f"{transcript_id}_exon{len(transcript_exons[transcript_id])+1}"
            
            # Extract exon_number
            exon_number_match = exon_number_pattern.search(attributes)
            exon_number = int(exon_number_match.group(1)) if exon_number_match else len(transcript_exons[transcript_id]) + 1
            
            # Store exon with its ENSE ID and exon_number
            transcript_exons[transcript_id].append((start, end, exon_id, exon_number))
            
            # Store strand and gene (only once per transcript)
            if transcript_id not in transcript_strand:
                transcript_strand[transcript_id] = strand
                transcript_gene[transcript_id] = gene_id
    
    if verbose:
        print(f"    Processed {line_count:,} lines total")
        print("  Sorting exons by exon_number...")
    
    # Sort exons by exon_number (which represents 5'->3' transcript order)
    for transcript_id in transcript_exons:
        transcript_exons[transcript_id].sort(key=lambda x: x[3])  # Sort by exon_number
    
    return transcript_exons, transcript_strand, transcript_gene


def extract_sequences(args, transcript_sequences: Dict[str, str], 
                    transcript_exons: Dict[str, List], 
                    transcript_strand: Dict[str, str],
                    transcript_gene: Dict[str, str]) -> bool:
    """Extract exon and junction sequences and write to output file."""
    
    # Define inadmissible junction configurations
    inadmissible_configs = [
        (22, 1),   # -22:1
        (21, 2),   # -21:2
        (20, 3),   # -20:3
        (19, 4),
        (18, 5),
        (17, 6),
        (6, 17),
        (5, 18),
        (4, 19),
        (3, 20),   # -3:20
        (2, 21),   # -2:21
        (1, 22),   # -1:22
    ]
    
    with open(args.output, 'w') as out:
        # Write header
        out.write("exon_id\ttranscript_id\tgene_id\texon_number\tstrand\texon_sequence\texon_length\tis_admissible\n")
        
        exon_count = 0
        junction_count = 0
        three_exon_count = 0
        transcripts_processed = 0
        transcripts_skipped = 0
        transcripts_skipped_length_mismatch = 0
        mismatch_count = 0
        
        for transcript_id, exons in transcript_exons.items():
            transcripts_processed += 1
            if args.verbose and transcripts_processed % 1000 == 0:
                print(f"  Processed {transcripts_processed:,} transcripts, {exon_count:,} exons, {junction_count:,} junctions...", end='\r')
            
            # Get transcript sequence
            if transcript_id not in transcript_sequences:
                transcripts_skipped += 1
                continue
            
            transcript_seq = transcript_sequences[transcript_id]
            gene_id = transcript_gene.get(transcript_id, "NA")
            strand = transcript_strand.get(transcript_id, "+")
            
            # VALIDATION: Check if GTF exon lengths sum to FASTA transcript length
            gtf_total_length = sum(exon_end - exon_start + 1 for exon_start, exon_end, _, _ in exons)
            fasta_length = len(transcript_seq)
            
            if gtf_total_length != fasta_length:
                if transcripts_skipped_length_mismatch < args.max_warnings:
                    print(f"\nWARNING: {transcript_id} length mismatch!")
                    print(f"  GTF exons sum to: {gtf_total_length:,} bp")
                    print(f"  FASTA sequence:   {fasta_length:,} bp")
                    print(f"  Difference:       {abs(gtf_total_length - fasta_length):,} bp")
                    print(f"  Skipping this transcript...")
                transcripts_skipped_length_mismatch += 1
                continue
            
            # Calculate cumulative positions in transcript sequence
            cumulative_pos = 0
            exon_positions = []  # Store (exon_id, start_pos, end_pos, exon_number) for junction extraction
            exon_lengths = []  # Track exon lengths to identify short exons
            
            for exon_start, exon_end, exon_id, exon_number in exons:
                # Calculate expected exon length from genomic coordinates
                expected_length = exon_end - exon_start + 1
                
                # Extract sequence from transcript
                exon_seq = transcript_seq[cumulative_pos:cumulative_pos + expected_length]
                actual_length = len(exon_seq)
                
                if actual_length != expected_length:
                    if mismatch_count < args.max_warnings:
                        print(f"\nWARNING: {transcript_id} exon {exon_number}: GTF says {expected_length}bp but extracted {actual_length}bp")
                    mismatch_count += 1
                    exon_length = actual_length
                else:
                    exon_length = expected_length
                
                # Write exon to output (unless junctions-only mode)
                if not args.junctions_only:
                    out.write(f"{exon_id}\t{transcript_id}\t{gene_id}\t{exon_number}\t{strand}\t{exon_seq}\t{exon_length}\tTrue\n")
                    exon_count += 1
                
                # Store exon position and length for junction extraction
                exon_positions.append((exon_id, cumulative_pos, cumulative_pos + exon_length, exon_number))
                exon_lengths.append(exon_length)
                
                cumulative_pos += exon_length
            
            # Skip junction processing if exons-only mode
            if args.exons_only:
                continue
            
            # Track which junction pairs have been processed for short exons
            # to avoid duplicate processing in the regular 2-exon loop
            short_exon_junctions: Set[Tuple[int, int]] = set()
            
            # First: Create 3-exon sequences for short middle exons
            for i in range(len(exon_positions) - 2):
                exon1_id, exon1_start, exon1_end, exon1_num = exon_positions[i]
                exon2_id, exon2_start, exon2_end, exon2_num = exon_positions[i + 1]
                exon3_id, exon3_start, exon3_end, exon3_num = exon_positions[i + 2]
                
                exon1_length = exon_lengths[i]
                exon2_length = exon_lengths[i + 1]
                exon3_length = exon_lengths[i + 2]
                
                # If middle exon (exon2) is short, create 3-exon sequence
                if exon2_length < args.short_exon_threshold:
                    # Track these junctions to handle specially
                    short_exon_junctions.add((i, i + 1))
                    short_exon_junctions.add((i + 1, i + 2))
                    
                    # Extract: last 19bp of exon1 + full exon2 + first 19bp of exon3
                    three_exon_seq = ""
                    
                    # Last 19bp of exon1 (or all if exon1 < 19bp)
                    if exon1_length >= 19:
                        three_exon_seq += transcript_seq[exon1_end - 19:exon1_end]
                    else:
                        three_exon_seq += transcript_seq[exon1_start:exon1_end]
                    
                    # Full exon2 (the short middle exon)
                    three_exon_seq += transcript_seq[exon2_start:exon2_end]
                    
                    # First 19bp of exon3 (or all if exon3 < 19bp)
                    if exon3_length >= 19:
                        three_exon_seq += transcript_seq[exon3_start:exon3_start + 19]
                    else:
                        three_exon_seq += transcript_seq[exon3_start:exon3_end]
                    
                    # Create 3-exon ID
                    three_exon_id = f"{exon1_id}-{exon2_id}-{exon3_id}"
                    three_exon_number = f"{exon1_num}-{exon2_num}-{exon3_num}"
                    
                    # Write 3-exon sequence (admissible)
                    out.write(f"{three_exon_id}\t{transcript_id}\t{gene_id}\t{three_exon_number}\t{strand}\t{three_exon_seq}\t{len(three_exon_seq)}\tTrue\n")
                    junction_count += 1
                    three_exon_count += 1
                    
                    # Create inadmissible 2-exon junctions with the short exon
                    if not args.skip_inadmissible:
                        # Junction: exon1 - exon2 (short)
                        junction_12_number = f"{exon1_num}-{exon2_num}"
                        for left_bp, right_bp in inadmissible_configs:
                            seq_12 = ""
                            if exon1_length >= left_bp:
                                seq_12 += transcript_seq[exon1_end - left_bp:exon1_end]
                            else:
                                seq_12 += transcript_seq[exon1_start:exon1_end]
                            
                            if exon2_length >= right_bp:
                                seq_12 += transcript_seq[exon2_start:exon2_start + right_bp]
                            else:
                                seq_12 += transcript_seq[exon2_start:exon2_end]
                            
                            junction_id_12 = f"{exon1_id}-{exon2_id}_inadmissible_{left_bp}:{right_bp}"
                            out.write(f"{junction_id_12}\t{transcript_id}\t{gene_id}\t{junction_12_number}\t{strand}\t{seq_12}\t{len(seq_12)}\tFalse\n")
                            junction_count += 1
                        
                        # Junction: exon2 (short) - exon3
                        junction_23_number = f"{exon2_num}-{exon3_num}"
                        for left_bp, right_bp in inadmissible_configs:
                            seq_23 = ""
                            if exon2_length >= left_bp:
                                seq_23 += transcript_seq[exon2_end - left_bp:exon2_end]
                            else:
                                seq_23 += transcript_seq[exon2_start:exon2_end]
                            
                            if exon3_length >= right_bp:
                                seq_23 += transcript_seq[exon3_start:exon3_start + right_bp]
                            else:
                                seq_23 += transcript_seq[exon3_start:exon3_end]
                            
                            junction_id_23 = f"{exon2_id}-{exon3_id}_inadmissible_{left_bp}:{right_bp}"
                            out.write(f"{junction_id_23}\t{transcript_id}\t{gene_id}\t{junction_23_number}\t{strand}\t{seq_23}\t{len(seq_23)}\tFalse\n")
                            junction_count += 1
            
            # Second: Extract regular 2-exon junction sequences
            for i in range(len(exon_positions) - 1):
                # Skip if this junction involves a short middle exon (already processed)
                if (i, i + 1) in short_exon_junctions:
                    continue
                
                exon1_id, exon1_start, exon1_end, exon1_num = exon_positions[i]
                exon2_id, exon2_start, exon2_end, exon2_num = exon_positions[i + 1]
                junction_exon_number = f"{exon1_num}-{exon2_num}"
                
                exon1_length = exon1_end - exon1_start
                exon2_length = exon2_end - exon2_start
                
                # TYPE 1: Admissible junction (-19:+19)
                junction_seq_admissible = ""
                
                # Last 19bp of exon1 (or all of it if shorter than 19bp)
                if exon1_length >= 19:
                    junction_seq_admissible += transcript_seq[exon1_end - 19:exon1_end]
                else:
                    junction_seq_admissible += transcript_seq[exon1_start:exon1_end]
                
                # First 19bp of exon2 (or all of it if shorter than 19bp)
                if exon2_length >= 19:
                    junction_seq_admissible += transcript_seq[exon2_start:exon2_start + 19]
                else:
                    junction_seq_admissible += transcript_seq[exon2_start:exon2_end]
                
                junction_exon_id_admissible = f"{exon1_id}-{exon2_id}"
                
                # Write admissible junction
                out.write(f"{junction_exon_id_admissible}\t{transcript_id}\t{gene_id}\t{junction_exon_number}\t{strand}\t{junction_seq_admissible}\t{len(junction_seq_admissible)}\tTrue\n")
                junction_count += 1
                
                # TYPES 2-7: Inadmissible junctions (asymmetric)
                if not args.skip_inadmissible:
                    for left_bp, right_bp in inadmissible_configs:
                        junction_seq_inadmissible = ""
                        
                        # Last N bp of exon1 (or all of it if shorter)
                        if exon1_length >= left_bp:
                            junction_seq_inadmissible += transcript_seq[exon1_end - left_bp:exon1_end]
                        else:
                            junction_seq_inadmissible += transcript_seq[exon1_start:exon1_end]
                        
                        # First M bp of exon2 (or all of it if shorter)
                        if exon2_length >= right_bp:
                            junction_seq_inadmissible += transcript_seq[exon2_start:exon2_start + right_bp]
                        else:
                            junction_seq_inadmissible += transcript_seq[exon2_start:exon2_end]
                        
                        # Include config in ID to make each inadmissible junction unique
                        junction_exon_id_inadmissible = f"{exon1_id}-{exon2_id}_inadmissible_{left_bp}:{right_bp}"
                        
                        # Write inadmissible junction
                        out.write(f"{junction_exon_id_inadmissible}\t{transcript_id}\t{gene_id}\t{junction_exon_number}\t{strand}\t{junction_seq_inadmissible}\t{len(junction_seq_inadmissible)}\tFalse\n")
                        junction_count += 1
        
        if args.verbose:
            print(f"\n  Processed {transcripts_processed:,} transcripts, {exon_count:,} exons, {junction_count:,} junctions total")
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"Exon and junction sequences written to: {args.output}")
    print(f"Total transcripts processed: {transcripts_processed:,}")
    print(f"Transcripts skipped (no FASTA): {transcripts_skipped:,}")
    print(f"Transcripts skipped (length mismatch): {transcripts_skipped_length_mismatch:,}")
    if not args.junctions_only:
        print(f"Total exons extracted: {exon_count:,}")
    if not args.exons_only:
        print(f"Total junctions extracted: {junction_count:,}")
        print(f"  - 3-exon sequences (short middle exons): {three_exon_count:,}")
    if mismatch_count > 0:
        print(f"Individual exon length mismatches: {mismatch_count:,}")
    print(f"{'='*60}")
    
    return True


def main():
    """Main entry point."""
    args = parse_args()
    
    print("Loading transcript sequences from FASTA...")
    transcript_sequences = load_fasta_sequences(args.fasta, args.verbose)
    print(f"  Loaded {len(transcript_sequences):,} transcript sequences")
    
    print("\nParsing GTF file for exons...")
    transcript_exons, transcript_strand, transcript_gene = parse_gtf_exons(args.gtf, args.verbose)
    print(f"  Found {len(transcript_exons):,} transcripts with exons")
    
    print("\nExtracting sequences...")
    success = extract_sequences(args, transcript_sequences, transcript_exons, 
                                transcript_strand, transcript_gene)
    
    return success


if __name__ == "__main__":
    start_time = time.time()
    
    try:
        success = main()
        elapsed = time.time() - start_time
        print(f"\nTotal runtime: {elapsed:.2f} seconds")
        
        if success:
            print("All done!")
            sys.exit(0)
        else:
            print("Error occurred")
            sys.exit(1)
    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
        sys.exit(130)
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)