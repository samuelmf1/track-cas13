#!/usr/bin/awk -f

BEGIN {
    FS = ","; OFS = ",";
}

NR == 1 {
    print $0                # header stdout
    print $0 > fails        # header fails
    next
}

# Good rows: third column does NOT contain the bad motifs
NR > 1 && $3 !~ /TTTTT|AAAA|CCCCC|GGGGG/ {
    print $0
    next
}

# Bad rows: redirect to fails
{
    print $0 > fails
}
