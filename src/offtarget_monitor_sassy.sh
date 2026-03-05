#!/bin/bash

BASE_DIR=$1
if [ -z "$BASE_DIR" ]; then
    echo "Usage: $0 <nextflow_project_root>"
    exit 1
fi

BASE_DIR=$(realpath "$BASE_DIR")
LOG_FILE="$BASE_DIR/.nextflow.log"
WORK_BASE="$BASE_DIR/work"
CACHE_FILE="/tmp/sassy_monitor_$(echo $BASE_DIR | md5sum | cut -d' ' -f1).cache"

if [ ! -f "$LOG_FILE" ]; then
    echo "Error: .nextflow.log not found at $LOG_FILE"
    exit 1
fi

# Header - Adjusted widths to accommodate (X/Y)
printf "%-5s %-20s %-10s %-25s\n" "ID" "Progress (%)" "Change" "Current Sequence"
printf "%-5s %-20s %-10s %-25s\n" "--" "------------" "------" "----------------"

CURRENT_RUN_CACHE=$(mktemp)

tac "$LOG_FILE" | grep "Submitted process > OFFTARGET:SASSY_SEARCH_CHUNK" | \
awk '{
    match($0, /\(([0-9]+)\)/, id); 
    match($0, /\[([a-z0-9/]+)\]/, wd);
    if (!seen[id[1]]++) print id[1] "|" wd[1]
}' | \
sort -n -t'|' -k1,1 | while IFS='|' read -r CHUNK_ID SHORT_WDIR; do

    WDIR="$WORK_BASE/$SHORT_WDIR"
    [[ ! -d "$WDIR" ]] && WDIR=$(find "$WORK_BASE" -maxdepth 2 -name "$(basename $SHORT_WDIR)*" -type d -print -quit)

    PERCENT="0.00"
    CURR_LABEL="Queued/Pending"
    COUNTS="0/0"

    if [ -d "$WDIR" ]; then
        FASTA_FILE="$WDIR/gencode.v49.filtered.part_${CHUNK_ID}.fa"
        OUTPUT_FILE="$WDIR/chunk_output.txt"

        if [[ -f "$FASTA_FILE" && -f "$OUTPUT_FILE" ]]; then
            LAST_SEQ=$(tail -n 1 "$OUTPUT_FILE" | awk '{print $2}')
            if [ ! -z "$LAST_SEQ" ]; then
                CURR_LINE=$(grep ">" "$FASTA_FILE" | grep -n -m 1 "$LAST_SEQ" | cut -d: -f1)
                TOTAL_LINES=$(grep -c ">" "$FASTA_FILE")
                
                if [[ ! -z "$CURR_LINE" && "$TOTAL_LINES" -gt 0 ]]; then
                    PERCENT=$(echo "scale=2; ($CURR_LINE / $TOTAL_LINES) * 100" | bc)
                    COUNTS="$CURR_LINE/$TOTAL_LINES"
                fi
                CURR_LABEL=$(echo "$LAST_SEQ" | cut -d'|' -f1)
            fi
        fi
    fi

    PREV_PERCENT=$(grep "^${CHUNK_ID}:" "$CACHE_FILE" 2>/dev/null | cut -d: -f2)
    DIFF="0.00"
    if [[ ! -z "$PREV_PERCENT" ]]; then
        DIFF=$(echo "$PERCENT - $PREV_PERCENT" | bc)
    fi

    if (( $(echo "$DIFF > 0" | bc -l) )); then
        DIFF_STR="+$DIFF%"
    else
        DIFF_STR="--"
    fi

    echo "${CHUNK_ID}:${PERCENT}" >> "$CURRENT_RUN_CACHE"
    
    # Combine Percent and Counts into one string
    PROGRESS_STR="${PERCENT}% (${COUNTS})"
    
    printf "%-5s %-20s %-10s %-25s\n" "$CHUNK_ID" "$PROGRESS_STR" "$DIFF_STR" "$CURR_LABEL"

done

mv "$CURRENT_RUN_CACHE" "$CACHE_FILE"