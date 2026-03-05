#!/usr/bin/awk -f

BEGIN {
    FS = ","
}

# Skip header row and process data
NR > 1 {
    printf ">%s|%s|%s|%s\n%s\n", $1, $2, $6, $7, $3
}