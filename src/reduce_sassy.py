#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(description="Filter guides based on Sassy blacklist")
    parser.add_argument("-i", "--input", required=True, help="Input guide CSV/TSV")
    parser.add_argument("-b", "--blacklist", required=True, help="Sassy blacklist (sequences only)")
    parser.add_argument("-o", "--output", required=True, help="Output filtered file")
    args = parser.parse_args()

    # Load input
    try:
        if args.input.endswith('.tsv'):
            df = pd.read_csv(args.input, sep='\t')
        else:
            df = pd.read_csv(args.input)
    except Exception as e:
        print(f"Error reading input {args.input}: {e}")
        sys.exit(1)

    # Load blacklist
    try:
        with open(args.blacklist, 'r') as f:
            blacklist = set(line.strip() for line in f if line.strip())
    except Exception as e:
        print(f"Error reading blacklist {args.blacklist}: {e}")
        sys.exit(1)

    if not blacklist:
        df.to_csv(args.output, index=False, sep='\t' if args.output.endswith('.tsv') else ',')
        return

    # Filter
    # In collapse.py output, the column is "Guide Sequence"
    if "Guide Sequence" not in df.columns:
        print(f"Warning: 'Guide Sequence' column not found in {args.input}. Skipping filtering.")
        df.to_csv(args.output, index=False, sep='\t' if args.output.endswith('.tsv') else ',')
        return

    filtered_df = df[~df["Guide Sequence"].isin(blacklist)]
    
    # Save output
    filtered_df.to_csv(args.output, index=False, sep='\t' if args.output.endswith('.tsv') else ',')
    print(f"Filtered {len(df) - len(filtered_df)} guides from {args.input}")

if __name__ == "__main__":
    main()
