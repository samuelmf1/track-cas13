#!/bin/bash
#SBATCH --job-name=sassy_orchestration
#SBATCH --output=/gpfs/commons/home/sfriedman/projects/tracklib/runtime/logs/sassy_%j_orchestration.out
#SBATCH --error=/gpfs/commons/home/sfriedman/projects/tracklib/runtime/logs/sassy_%j_orchestration.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=48:00:00

INPUT_FILE=$1
NUM_CHUNKS=7500  # You can change this to 100, 500, etc.
BASENAME=$(basename "$INPUT_FILE")
TMP_DIR="tmp_${BASENAME}_$(date +%s)"
mkdir -p "$TMP_DIR"

# 1. Calculate lines per chunk
TOTAL_LINES=$(wc -l < "$INPUT_FILE")

# If the file is smaller than our target chunk count, 
# we just do 1 line per file.
if [ "$TOTAL_LINES" -le "$NUM_CHUNKS" ]; then
    LPP=1
else
    # Calculate lines per partition, rounding up
    LPP=$(( (TOTAL_LINES + NUM_CHUNKS - 1) / NUM_CHUNKS ))
fi

# 2. Split the file
# -d: numeric suffixes, -a 4: allows up to 9999 chunks
split -l "$LPP" -d -a 4 "$INPUT_FILE" "$TMP_DIR/chunk_"

# 3. Determine the actual upper bound of the array
# We count how many files were actually created
ACTUAL_CHUNKS=$(ls "$TMP_DIR"/chunk_* | wc -l)
MAX_INDEX=$((ACTUAL_CHUNKS - 1))

echo "Split $TOTAL_LINES lines into $ACTUAL_CHUNKS chunks ($LPP lines each)."

# 4. Submit the Array Job
# Passing TMP_DIR so the worker knows where to look
ARRAY_JOB_ID=$(sbatch --parsable --array=0-$MAX_INDEX --export=ALL,TMP_DIR="$TMP_DIR" /gpfs/commons/home/sfriedman/projects/tracklib/src/sassy_worker.sh)

# 5. Submit the Merge Job (Wait for the array to finish)
sbatch --dependency=afterok:$ARRAY_JOB_ID --export=ALL,TMP_DIR="$TMP_DIR",OUT_NAME="$BASENAME" /gpfs/commons/home/sfriedman/projects/tracklib/src/sassy_merge.sh

echo "Orchestration started. Array ID: $ARRAY_JOB_ID"