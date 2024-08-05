#!/usr/bin/env python3
# coding=utf-8

import argparse
import datetime
from pathlib import Path
import pandas as pd

# TODO add processing of failed and putative files
# TODO add merge with previous retrofit 

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
    
    merged_df['Action'] = 'add'

    result_df = merged_df[['Primary_assembly', 'MAG_accession', 'Genome', 'Species_rep', 'Action']]

    if previous_table:
        previous_df = pd.read_csv(previous_table, sep='\t')
        unique_previous_df = previous_df[~previous_df['MAG_accession'].isin(result_df['MAG_accession'])]
        result_df = pd.concat([unique_previous_df, result_df]).reset_index(drop=True)
        
    result_df.fillna('NA', inplace=True)
    output_file = generate_filename("mag_to_assembly_links")
    result_df.to_csv(output_file, sep='\t', index=False)


def generate_filename(prefix):
    current_date = datetime.datetime.now().strftime(("%Y-%m-%d_%Hh%Mm"))
    return f"{prefix}_{current_date}.tsv"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create assembly-genome linking table from the output files of the main script")
    parser.add_argument("--previous-table", 
                        "-p",
                        default=None, 
                        help="Linking table generated in the previous run of the pipeline")
    parser.add_argument("--catalogue-metadata", 
                        "-m" ,
                        required=True,
                        type=Path, 
                        help="Metadata file from MGnify catalogues to take Species_rep IDs")
    parser.add_argument('results', 
                        metavar='FILE', 
                        nargs='+',
                        help="Files with links produced by the main script")
    args = parser.parse_args()
    main(args.results, args.catalogue_metadata, args.previous_table)