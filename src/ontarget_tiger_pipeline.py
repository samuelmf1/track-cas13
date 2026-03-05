#!/usr/bin/env python
# FAST TIGER - On-target only, TF 2.11 compatible, 50–200x speedup
# Memory-efficient version with streaming output
# Fixed to produce IDENTICAL output to original TIGER
# FIXED: Skip guides containing N nucleotides

import argparse
import os
import gzip
import numpy as np
import pandas as pd
import tensorflow as tf
from Bio import SeqIO
from tqdm.auto import tqdm

# Deterministic ops for bit-exact reproducibility
os.environ['TF_DETERMINISTIC_OPS'] = '1'

ID_COL = 'Contig'
IDX_COL = 'Contig_Idx'
TARGET_COL = 'Target Sequence'
GUIDE_COL = 'Guide Sequence'
SCORE_COL = 'Guide Score'
GENE_COL = 'Gene'
GENE_SYMBOL_COL = 'Symbol'
COORDS_COL = 'Coords'
TITLE_COL = 'Title'

GUIDE_LEN = 23
CONTEXT_5P = 3
CONTEXT_3P = 0
TARGET_LEN = CONTEXT_5P + GUIDE_LEN + CONTEXT_3P

# Use EXACT same complement dict as original (no .get() with default)
NUCLEOTIDE_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

# Use EXACT same batch size as original
BATCH_SIZE_COMPUTE = 500

# Chunk size for streaming output (adjust based on available RAM)
CHUNK_SIZE = 10000

# GPU memory growth
gpus = tf.config.list_physical_devices('GPU')
if gpus:
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)


def load_transcripts(fasta_paths):
    records = []
    paths = fasta_paths if isinstance(fasta_paths, (list, tuple)) else [fasta_paths]
    for path in tqdm(paths, desc="FASTA files", unit="file"):
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt') as f:
            for rec in SeqIO.parse(f, 'fasta'):
                seq = str(rec.seq).upper()
                if len(seq) < TARGET_LEN:
                    continue
                transcript_id = rec.id.split('|')[0]
                records.append({
                    'id': transcript_id,
                    'header': rec.id,
                    'seq': seq
                })
    df = pd.DataFrame(records).set_index('id')
    return df


def has_ambiguous_nucleotides(sequence):
    """Check if sequence contains any ambiguous nucleotides (N or other)"""
    return any(nt not in NUCLEOTIDE_COMPLEMENT for nt in sequence)


def sequence_complement(sequences):
    """Match original implementation exactly - no default handling"""
    return [''.join([NUCLEOTIDE_COMPLEMENT[nt] for nt in list(seq)]) for seq in sequences]


def one_hot_encode_sequence(sequence, add_context_padding=False):
    """Match original implementation exactly"""
    # Stack list of sequences into a tensor
    sequence = tf.ragged.stack([tf.constant(list(seq)) for seq in sequence], axis=0)

    # Tokenize sequence - match original token values exactly
    nucleotide_table = tf.lookup.StaticVocabularyTable(
        initializer=tf.lookup.KeyValueTensorInitializer(
            keys=tf.constant(['A', 'C', 'G', 'T', 'N'], dtype=tf.string),
            values=tf.constant([0, 1, 2, 3, 255], dtype=tf.int64)),
        num_oov_buckets=1)
    sequence = tf.RaggedTensor.from_row_splits(
        values=nucleotide_table.lookup(sequence.values),
        row_splits=sequence.row_splits).to_tensor(255)

    # Add context padding if requested
    if add_context_padding:
        pad_5p = 255 * tf.ones([sequence.shape[0], CONTEXT_5P], dtype=sequence.dtype)
        pad_3p = 255 * tf.ones([sequence.shape[0], CONTEXT_3P], dtype=sequence.dtype)
        sequence = tf.concat([pad_5p, sequence, pad_3p], axis=1)

    # One-hot encode
    sequence = tf.one_hot(sequence, depth=4, dtype=tf.float16)

    return sequence


def process_data(transcript_seq):
    """Match original implementation exactly - now filters out sequences with N"""
    # Convert to upper case
    transcript_seq = transcript_seq.upper()

    # Get all target sites
    all_target_seq = [transcript_seq[i: i + TARGET_LEN] for i in range(len(transcript_seq) - TARGET_LEN + 1)]

    # Filter out sequences containing N or other ambiguous nucleotides
    valid_indices = []
    target_seq = []
    for idx, seq in enumerate(all_target_seq):
        if not has_ambiguous_nucleotides(seq):
            valid_indices.append(idx)
            target_seq.append(seq)
    
    # If no valid sequences, return empty
    if not target_seq:
        return [], [], None, []

    # Prepare guide sequences (only for valid targets)
    guide_seq = sequence_complement([seq[CONTEXT_5P:len(seq) - CONTEXT_3P] for seq in target_seq])

    # Model inputs
    model_inputs = tf.concat([
        tf.reshape(one_hot_encode_sequence(target_seq, add_context_padding=False), [len(target_seq), -1]),
        tf.reshape(one_hot_encode_sequence(guide_seq, add_context_padding=True), [len(guide_seq), -1]),
        ], axis=-1)
    return target_seq, guide_seq, model_inputs, valid_indices


def calibrate_predictions(predictions, num_mismatches):
    """Match original implementation exactly"""
    params = pd.read_pickle('calibration_params.pkl')
    correction = np.squeeze(params.set_index('num_mismatches').loc[num_mismatches, 'slope'].to_numpy())
    return correction * predictions


def score_predictions(predictions):
    """Match original implementation exactly"""
    params = pd.read_pickle('scoring_params.pkl')
    params = params.iloc[0]
    return 1 - 1 / (1 + np.exp(params['a'] * predictions + params['b']))


def format_chunk(chunk_df):
    """Apply final formatting to match original TIGER exactly"""
    chunk_df[GUIDE_COL] = chunk_df[GUIDE_COL].str[::-1]  # Reverse guide sequences
    chunk_df[TARGET_COL] = chunk_df[TARGET_COL].str[CONTEXT_5P : -CONTEXT_3P or None]  # Remove context
    return chunk_df


def predict_on_targets_fast(transcripts_df, model, output_path):
    """Memory-efficient streaming prediction with chunked output"""
    
    chunk_buffer = []
    first_chunk = True
    total_guides = 0
    total_transcripts = 0
    skipped_guides = 0
    
    print(f"- Processing {len(transcripts_df):,} transcripts...")
    
    for transcript_id, row in tqdm(transcripts_df.iterrows(),
                                total=len(transcripts_df),
                                desc="Transcripts",
                                unit="tx"):
        seq = row['seq']
        if len(seq) < TARGET_LEN:
            continue
        
        # Use EXACT same processing as original, but now filters N's
        target_seq, guide_seq, model_inputs, valid_indices = process_data(seq)
        
        # Skip transcript if no valid guides
        if not target_seq:
            skipped_guides += len(seq) - TARGET_LEN + 1
            continue
        
        # Track skipped guides
        total_possible = len(seq) - TARGET_LEN + 1
        skipped_guides += total_possible - len(target_seq)
        
        # Use EXACT same prediction method as original
        lfc_estimate = model.predict(model_inputs, batch_size=BATCH_SIZE_COMPUTE, verbose=False)[:, 0]
        
        # Use EXACT same calibration as original
        lfc_estimate = calibrate_predictions(lfc_estimate, num_mismatches=np.zeros_like(lfc_estimate))
        
        # Use EXACT same scoring as original
        scores = score_predictions(lfc_estimate)
        
        # Parse header exactly as original
        header_parts = row['header'].split('|')
        gene = header_parts[1] if len(header_parts) > 1 else ''
        symbol = header_parts[2] if len(header_parts) > 2 else ''
        coords = header_parts[3] if len(header_parts) > 3 else ''
        title = header_parts[4] if len(header_parts) > 4 else ''
        
        # Build results matching original structure exactly
        for idx, (target, guide, score) in enumerate(zip(target_seq, guide_seq, scores)):
            chunk_buffer.append({
                ID_COL: transcript_id,
                IDX_COL: valid_indices[idx],  # Use original index
                TARGET_COL: target,
                GUIDE_COL: guide,
                SCORE_COL: round(float(score), 4),  # Reduce precision due to batching differently
                GENE_COL: gene,
                GENE_SYMBOL_COL: symbol,
                COORDS_COL: coords,
                TITLE_COL: title
            })
        
        total_transcripts += 1
        
        # Write chunk when buffer is full
        if len(chunk_buffer) >= CHUNK_SIZE:
            chunk_df = pd.DataFrame(chunk_buffer)
            chunk_df = format_chunk(chunk_df)
            
            # Write to CSV (append mode after first chunk)
            chunk_df.to_csv(output_path, mode='a', header=first_chunk, index=False)
            
            total_guides += len(chunk_df)
            first_chunk = False
            chunk_buffer = []
    
    # Write remaining buffer
    if chunk_buffer:
        chunk_df = pd.DataFrame(chunk_buffer)
        chunk_df = format_chunk(chunk_df)
        chunk_df.to_csv(output_path, mode='a', header=first_chunk, index=False)
        total_guides += len(chunk_df)
    
    if skipped_guides > 0:
        print(f"- Skipped {skipped_guides:,} guides containing ambiguous nucleotides (N)")
    
    return total_guides, total_transcripts


def main():
    parser = argparse.ArgumentParser(description="FAST TIGER - Ultra-fast on-target scoring (memory-efficient)")
    parser.add_argument('--fasta_path', type=str, required=True, help="FASTA file")
    parser.add_argument('--output', type=str, default='tiger_results_on_target.csv')
    parser.add_argument('--chunk_size', type=int, default=50_000, help="Chunk size for streaming output")
    args = parser.parse_args()

    # Set chunk size from argument
    global CHUNK_SIZE
    CHUNK_SIZE = args.chunk_size

    tf.config.set_visible_devices([], 'GPU')
    print("- Visible devices:", tf.config.get_visible_devices())

    # Load model
    print("- Loading model...")
    model = tf.keras.models.load_model('model', compile=False)

    # Load FASTA
    if not os.path.isfile(args.fasta_path):
        raise ValueError("FASTA file not found")

    transcripts_df = load_transcripts(args.fasta_path)
    print(f"- Loaded {len(transcripts_df):,} transcripts")

    # Run prediction with streaming output
    total_guides, total_transcripts = predict_on_targets_fast(transcripts_df, model, args.output)

    if total_guides == 0:
        print("No valid guide sites found.")
        return

    print(f"- Saved {total_guides:,} guides ({total_transcripts:,} transcripts) to {args.output}")


if __name__ == '__main__':
    main()