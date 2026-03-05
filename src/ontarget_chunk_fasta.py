#!/usr/bin/env python3
from Bio import SeqIO
import math
import os
import argparse

def split_fasta(input_fasta, start_seq_id=1, n_chunks=4):
    # Get directory of input file
    input_dir = os.path.dirname(input_fasta) or "."  # "." if input_fasta has no path
    input_basename = os.path.basename(input_fasta).replace(".fasta", "")
    
    # # Output prefix: {input_dir}/chunks/{basename}_chunked
    # chunks_dir = os.path.join(input_dir, "chunks")
    # output_fasta_prefix = f"{input_basename}_chunked" #os.path.join(chunks_dir, f"{input_basename}_chunked")
    
    # Collect sequences starting from start_seq_id
    records_to_write = [
        record for record in SeqIO.parse(input_fasta, "fasta")
        if int(record.id.split("|")[0]) >= start_seq_id
    ]
    
    total_records = len(records_to_write)
    if total_records == 0:
        print("No sequences found with the given start_seq_id.")
        return

    # Group records by gene_id (index 1 in pipe-delimited header)
    gene_groups = {} # gene_id -> list of records
    # We want to maintain order of appearance for gene_ids
    ordered_gene_ids = []
    
    for record in records_to_write:
        parts = record.id.split("|")
        # Ensure we have enough parts; if not, fallback or error? 
        # Assuming format is consistent as per user request example: >1|ENSG...|...
        if len(parts) > 1:
            gene_id = parts[1]
        else:
            gene_id = "unknown"
            
        if gene_id not in gene_groups:
            gene_groups[gene_id] = []
            ordered_gene_ids.append(gene_id)
        gene_groups[gene_id].append(record)

    # Calculate target chunk size (number of records per chunk)
    chunk_size_target = math.ceil(total_records / n_chunks)
    
    chunks = []
    current_chunk = []
    current_chunk_size = 0
    
    for gene_id in ordered_gene_ids:
        records = gene_groups[gene_id]
        
        # If adding this group exceeds target significantly, and we aren't at the last chunk,
        # we might want to split? BUT requirement is "all shared ENSG gene_ids all end up in the same chunk"
        # So we MUST add the whole group to the current chunk or the next one.
        
        # Simple greedy approach: fill current chunk until it hits target, then start next.
        # But if we are just starting a new chunk, we must add it.
        
        # If current chunk is already "full" (>= target) and we have more chunks to make, start new one.
        if current_chunk_size >= chunk_size_target and len(chunks) < n_chunks - 1:
            chunks.append(current_chunk)
            current_chunk = []
            current_chunk_size = 0
            
        current_chunk.extend(records)
        current_chunk_size += len(records)
        
    # Add the last chunk
    if current_chunk:
        chunks.append(current_chunk)

    # Write chunks
    for i, chunk_records in enumerate(chunks):
        if not chunk_records:
            continue
            
        output_file = f"{input_basename}.chunk_{i+1:02d}.fasta"
        SeqIO.write(chunk_records, output_file, "fasta")
        print(f"Wrote {len(chunk_records)} sequences to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Split a FASTA file into N chunks starting from a specific sequence ID."
    )
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument(
        "--start_seq_id", type=int, default=1,
        help="Numeric ID to start splitting from (default: 1)"
    )
    parser.add_argument(
        "--n_chunks", type=int, default=4,
        help="Number of chunks to split into (default: 4)"
    )

    args = parser.parse_args()
    split_fasta(args.input_fasta, args.start_seq_id, args.n_chunks)


if __name__ == "__main__":
    main()
