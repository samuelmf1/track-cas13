#!/bin/bash
#SBATCH --job-name=sassy_worker
#SBATCH --output=/gpfs/commons/home/sfriedman/projects/tracklib/runtime/logs/sassy_%j_worker.out
#SBATCH --error=/gpfs/commons/home/sfriedman/projects/tracklib/runtime/logs/sassy_%j_worker.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=3G
#SBATCH --time=24:00:00

FASTA_REF="/gpfs/commons/home/sfriedman/projects/tracklib/datasets/index/gencode.v38.filtered.fa"
SASSY_BIN="/gpfs/commons/home/sfriedman/.cargo/bin/sassy"

# Match the 4-digit suffix from the split command
CHUNK_ID=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
PATTERN_FILE="$TMP_DIR/chunk_$CHUNK_ID"

if [ -f "$PATTERN_FILE" ]; then
    $SASSY_BIN search \
        -j $SLURM_CPUS_PER_TASK \
        -k 2 --no-rc --max-n-frac 0 \
        --pattern-file "$PATTERN_FILE" \
        "$FASTA_REF" > "$TMP_DIR/result_$CHUNK_ID.txt"
else
    echo "Error: $PATTERN_FILE not found."
    exit 1
fi