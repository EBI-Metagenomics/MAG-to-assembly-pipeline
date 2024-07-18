#!/usr/bin/env python3
# coding=utf-8

import argparse
from pathlib import Path
import pandas as pd

def main(results, catalogue_metadata, previous_table):
    df_list = []
    for file in results:
        try:
            df = pd.read_csv(file, sep='\t', usecols=[0, 3], names=['MAG_accession', 'Primary_assembly'], header=None)
            df_list.append(df)
        except Exception as e:
            print(f"Error reading {file}: {e}")

    combined_df = pd.concat(df_list, ignore_index=True)
    expanded_df = combined_df.assign(Primary_assembly=combined_df['Primary_assembly'].str.split(',')).explode('Primary_assembly')

    metadata_df = pd.read_csv(catalogue_metadata, sep='\t')
    merged_df = pd.merge(expanded_df, metadata_df[['Genome', 'Species_rep', "Genome_accession"]],
                     left_on='MAG_accession', right_on='Genome_accession', how='left')
    merged_df.drop(columns=['Genome_accession'], inplace=True)

    merged_df.fillna('NA', inplace=True)

    merged_df['Action'] = 'add'

    final_df = merged_df[['Primary_assembly', 'MAG_accession', 'Genome', 'Species_rep', 'Action']]

    final_df.to_csv('RETROFIT.tsv', sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create assembly-genome linking table from the output files of the main script")
    parser.add_argument("--previous-table", 
                        "-p",
                        default=None, 
                        help="")
    parser.add_argument("--catalogue-metadata", 
                        "-m" ,
                        required=True,
                        type=Path, 
                        help="")
    parser.add_argument('results', 
                        metavar='FILE', 
                        nargs='+',
                        help="Files produced by the main script")
    args = parser.parse_args()
    main(args.results, args.catalogue_metadata, args.previous_table)