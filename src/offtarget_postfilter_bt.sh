#!/bin/bash
set -euo pipefail  # Exit on error, undefined vars, or pipe failures

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <full_data_file> <pairs_file> <output_prefix>" >&2
    exit 1
fi

FULL_DATA="$1"
PAIRS="$2"
OUTPUT_PREFIX="$3"

MATCHING="${OUTPUT_PREFIX}.bt.ot.csv"
NOT_MATCHING="${OUTPUT_PREFIX}.bt.csv"

# Use awk to split the data
# We use 'unbuffer' logic or just standard redirect
awk -F',' -v pairs_file="$PAIRS" -v matching="$MATCHING" -v not_matching="$NOT_MATCHING" '
BEGIN {
    while ((getline < pairs_file) > 0) {
        # Using the first two columns as the unique alignment key
        pairs[$1 "," $2] = 1
    }
    close(pairs_file)
}
NR == 1 {
    header = $0
    print header > matching
    print header > not_matching
    next
}
{
    key = $1 "," $2
    if (key in pairs) {
        print $0 > matching
    } else {
        print $0 > not_matching
    }
}
' "$FULL_DATA"

# Optional: Log stats to stderr so they appear in Nextflow's .command.err
echo "Stats for ${OUTPUT_PREFIX}: Total: $(wc -l < "$FULL_DATA"), Off-Target: $(wc -l < "$MATCHING")" >&2