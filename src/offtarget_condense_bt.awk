#!/usr/bin/awk -f
BEGIN { 
    FS = "\t"
    OFS = ","
    if (mode != "tx" && mode != "gene") {
        print "Error: mode must be 'tx' or 'gene'" > "/dev/stderr"
        exit 1
    }
    is_tx = (mode == "tx")
}

{
    split($1, q, "|")
    split($3, r, "|")
    
    if (is_tx ? (q[4] != r[1]) : (q[3] != r[2])) {
        key = q[1] SUBSEP q[2]
        if (!(key in seen)) {
            seen[key] = 1
            print q[1], q[2]
        }
    }
}