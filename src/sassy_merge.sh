#!/bin/bash
#SBATCH --job-name=sassy_merge
#SBATCH --output=/gpfs/commons/home/sfriedman/projects/tracklib/runtime/logs/sassy_%j_merge.out
#SBATCH --error=/gpfs/commons/home/sfriedman/projects/tracklib/runtime/logs/sassy_%j_merge.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

# Merge all results into the final blacklist
cat $TMP_DIR/result_*.txt > "${OUT_NAME}.blacklist.txt"

# Clean up temporary directory
rm -rf "$TMP_DIR"

echo "Processing complete: ${OUT_NAME}.blacklist.txt"