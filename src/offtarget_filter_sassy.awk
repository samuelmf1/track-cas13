#!/usr/bin/awk -f

BEGIN { 
    OFS="," 
}

# PHASE 1: Load the clean blacklist into memory
# The blacklist is the first file (NR==FNR)
# Expected format: [Sequence] [ENSG_ID]
NR == FNR {
    # $1 is sequence, $2 is ENSG ID
    bad_guides[$1] = $2
    next
}

# PHASE 2: Filter the Bowtie CSV (Comma-separated)
# We detect the first line of the second file to switch the delimiter
FNR == 1 { 
    FS=","
    $0=$0 # Re-parse line with new FS
    print $0
    next 
}

{
    current_gene = $1
    current_seq  = $3

    # If sequence is blacklisted...
    if (current_seq in bad_guides) {
        # ...AND the gene is NOT the target gene, drop it
        if (current_gene != bad_guides[current_seq]) {
            next
        }
    }
    
    print $0
}