#!/usr/bin/env python
import pandas as pd
import argparse
from tqdm import tqdm
try:
    from select_non_overlapping import GuideSelector, SelectionConfig
except ImportError:
    import sys
    import os
    # If running from src dir, the module name is select_constitutive
    sys.path.append(os.path.dirname(os.path.realpath(__file__)))
    try:
        from select_constitutive import GuideSelector, SelectionConfig
    except ImportError:
        # Final fallback: relative to current dir
        sys.path.append(os.getcwd())
        try:
            from select_non_overlapping import GuideSelector, SelectionConfig
        except ImportError:
            # If all else fails, pandas might be missing or name is different
            pass

def main():
    parser = argparse.ArgumentParser(description="Select guides for transcripts with unique targeting.")
    parser.add_argument("--tsv", type=str, required=True, help="Input TSV/CSV file containing guides")
    parser.add_argument("--output", type=str, required=True, help="Output TSV file")
    args = parser.parse_args()

    tqdm.pandas()

    # 1. Load Data
    df = pd.read_csv(args.tsv, sep=None, engine='python')

    # 2. Map columns (standardize)
    rename_map = {
        'Gene': 'gene_id',
        'Symbol': 'transcript_id_group',
        'Guide Score': 'tiger_score',
        'Contig_Idx': 'position',
        'Title': 'biotype'
    }
    df = df.rename(columns=rename_map)

    # 3. Ensure required columns exist
    for col in ['guide_expression_norm', 'ntargeted_tx', 'tags', 'biotype']:
        if col not in df.columns:
            if col == 'ntargeted_tx':
                df[col] = df['transcript_id_group'].apply(lambda x: str(x).count('|') + 1)
            elif col == 'tags':
                df[col] = ''
            elif col == 'biotype':
                df[col] = "NA"
            else:
                df[col] = 0

    # 4. FILTER: Keep only guides catering to a single transcript (Unique Targeting)
    df_unique = df[df['ntargeted_tx'] == 1].copy()
    
    # 5. Select Best Group AND Guides
    if df_unique.empty:
        select_df = pd.DataFrame(columns=df.columns)
    else:
        # Import here to catch issues if needed
        try:
            from select_non_overlapping import GuideSelector, SelectionConfig
            selector = GuideSelector(target_n=5)
            config = SelectionConfig(target_n=5)
            
            selected_dfs = []
            for gene_id, group in df_unique.groupby('gene_id'):
                res_df, _ = selector.process_gene(gene_id, group, config)
                if not res_df.empty:
                    selected_dfs.append(res_df)
            
            if selected_dfs:
                select_df = pd.concat(selected_dfs)
            else:
                select_df = pd.DataFrame(columns=df.columns)
        except ImportError as e:
            print(f"Error importing selection logic: {e}")
            # Fallback to empty if logic is missing
            select_df = pd.DataFrame(columns=df.columns)

    # 7. Final Output
    select_df.to_csv(args.output, sep="\t", index=False)
    print(f"Done. Selected {len(select_df)} guides.")

if __name__ == "__main__":
    main()
