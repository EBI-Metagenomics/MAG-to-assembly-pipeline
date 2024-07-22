#!/usr/bin/env python3
# coding=utf-8

import argparse
from ftplib import FTP
import logging 
import os
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        # logging.FileHandler(filename='script.log'),
        logging.StreamHandler()
    ]
)


mag_layer_file = Path("mag_layer.tsv")
bin_layer_file = Path("bin_layer.tsv")


def main(processed_acc_file, output_file, catalogues_metadata_file, gut_mapping_file):
    logging.info("Starting preparation of the list of input accessions for genome-primary assembly linking script...")
    logging.info("Step 1/5:")
    logging.info("Download all genomes from MGnify catalogues...")
    download_all_catalogues_metadata(Path("."), catalogues_metadata_file)
    remove_gut_genomes(gut_mapping_file, catalogues_metadata_file)

    logging.info("Step 2/5:")
    logging.info("Download all genomes from MAG layer of ENA...")
    mag_layer_data = {
        "result": "wgs_set",
        "query": 'assembly_type="metagenome-assembled genome (mag)"',
        "fields": """
                    wgs_set,
                    assembly_type,
                    assembly_accession,
                    sample_accession,
                    secondary_sample_accession,
                    run_accession,
                    study_accession,
                    status,
                    set_fasta_ftp
                    """,
        "format": "tsv"
    }
    download_layer(mag_layer_data, mag_layer_file)

    logging.info("Step 3/5:")
    logging.info("Download all genomes from bin layer of ENA...")
    bin_layer_data = {
        "result": "analysis",
        "query": 'assembly_type="binned metagenome"',
        "fields": """
                    analysis_accession,
                    assembly_type,
                    analysis_type,
                    sample_accession,
                    secondary_sample_accession,
                    run_accession,
                    secondary_study_accession,
                    study_accession,
                    status,
                    submitted_ftp,
                    generated_ftp
                    """ ,
        "format": "tsv"
    }
    download_layer(bin_layer_data, bin_layer_file)

    logging.info("Step 4/5:")
    logging.info("Merge all downloaded accessions removing redunduncy...")  
    catalogues_df = pd.read_csv(catalogues_metadata_file, usecols=[2], names=['Genome_accession'],  sep='\t')
    mag_layer_df = pd.read_csv(mag_layer_file, usecols=[0,2], names=['Genome_accession', 'Genome_NCBI_accession'], sep='\t')
    bin_layer_df = pd.read_csv(bin_layer_file, usecols=[0], names=['Genome_accession'], sep='\t')
    ncbi_accessions = catalogues_df[catalogues_df['Genome_accession'].str.startswith('GCA')]['Genome_accession']
    mag_layer_df = mag_layer_df[~mag_layer_df['Genome_NCBI_accession'].isin(ncbi_accessions)]
    combined_df = pd.concat([catalogues_df, mag_layer_df[['Genome_accession']], bin_layer_df], ignore_index=True)
    combined_df = combined_df.drop_duplicates()

    logging.info("Step 5/5:")
    logging.info("Remove from the output all accessions that found in the input list of previously processed genomes...")
    if processed_acc_file:
        processed_df = pd.read_csv(processed_acc_file, header=None, names=["Genome_accession"], sep='\t')
        combined_df = combined_df[~combined_df['Genome_accession'].isin(processed_df['Genome_accession'])]
    else: 
        logging.info("There is no list provided, skipping.")
    combined_df.to_csv(output_file, sep='\t', header=False, index=False)

    logging.info(f"Finished successfully! List of genome accessions is saved to {output_file}")


def download_all_catalogues_metadata(work_dir, output_file, previous=None):
    logging.info("Collecting all genomes-all_metadata.tsv files from MGnify FTP...")
    ftp_host = "ftp.ebi.ac.uk"
    ftp_dir = "/pub/databases/metagenomics/mgnify_genomes/"

    ftp = FTP(ftp_host)
    ftp.login()
    ftp.cwd(ftp_dir)
    catalogs = ftp.nlst()
    df_list = []
    for folder in tqdm(catalogs):
        ftp.cwd(folder)
        catalog_versions = ftp.nlst()
        latest_folder = which_is_latest(catalog_versions)
        if latest_folder:

            # if the catalog was processed previously, then skip it
            if previous and folder in previous and previous[folder] == latest_folder:
                continue

            ftp.cwd(latest_folder)
            filename = 'genomes-all_metadata.tsv'
            local_filename = work_dir / filename
            with open(local_filename, 'wb') as f:
                ftp.retrbinary('RETR ' + filename, f.write)
            df = pd.read_csv(local_filename, sep='\t')
            df = df[df["Genome_type"] == "MAG"]
            df = df[["Genome","Genome_type","Genome_accession","Species_rep"]]
            df_list.append(df)
            # TODO return this variable
            if previous:
                previous[folder] = latest_folder
            os.remove(local_filename)
            ftp.cwd('../..')
    ftp.quit()
    if df_list:
        concatenated_df = pd.concat(df_list, ignore_index=True)
        concatenated_df.to_csv(output_file.with_suffix('.source.tsv'), sep='\t', index=False)
        logging.info("Successfully completed!")
    else:
        logging.info("No catalogue updates were made since the last run of the pipeline. Nothing was downloaded.")


def which_is_latest(list_of_versions):
    latest_version = "0.0.0"
    latest_folder = None
    for version in list_of_versions:
        if version.startswith('v'):
            version_parts = version.rstrip('/').split('v')[-1].split('.')
            version_parts = [int(part) for part in version_parts]
            version_parts += [0] * (3 - len(version_parts))
            if tuple(version_parts) > tuple(map(int, latest_version.split('.'))):
                latest_version = '.'.join(map(str, version_parts))
                latest_folder = version
    return latest_folder


def remove_gut_genomes(gut_mapping_file, catalogues_metadata_file):
    logging.info("Starting preprocessing of MGnify catalogues metadata to remove GUT_GENOME ids...")
    catalogues_metadata_df = pd.read_csv(catalogues_metadata_file.with_suffix('.source.tsv'), sep='\t')
    gut_mapping_df = pd.read_csv(gut_mapping_file, sep='\t')
    replacement_dict = dict(zip(gut_mapping_df['Genome'], gut_mapping_df['Genome_accession']))

    def replace_genome_accession(row):
        if row['Genome_accession'].startswith('GUT_GENOME'):
            return replacement_dict.get(row['Genome'], None)
        return row['Genome_accession']

    catalogues_metadata_df['Genome_accession'] = catalogues_metadata_df.apply(replace_genome_accession, axis=1)
    catalogues_metadata_df = catalogues_metadata_df.dropna(subset=['Genome_accession'])
    catalogues_metadata_df.to_csv(catalogues_metadata_file, sep='\t', index=False)
    logging.info(f"Preprocessing completed successfully. Output saved to {catalogues_metadata_file}.")


def download_layer(data, output_file):
    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    try:
        response = requests.post(url, data=data, headers={"Content-Type": "application/x-www-form-urlencoded"})
        response.raise_for_status()  # Raise an error for bad status codes
        with open(output_file, 'w') as file:
            file.write(response.text)
        logging.info(f"Data downloaded to {output_file}")
    except requests.RequestException as e:
        logging.error(f"Failed to download data: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download lists of accessions from MGnify catalogues and ENA, merge them, remove redundancy and accessions that were processed in the past")
    parser.add_argument('--processed-acc', 
                        '-i',
                        default=None,
                        type=Path,
                        help='File with a list of accessions that were previously processed by the pipeline to skip them in the current run')
    parser.add_argument('--gut-mapping', 
                        '-g',
                        required=True,
                        type=Path,
                        help='File with a mapping of GUT_GENOME* accessions from human gut catalog too their source ENA accessions')
    parser.add_argument('--catalogue-metadata', 
                        '-m',
                        required=True,
                        type=Path,
                        help='File to write metadata of used accessions from MGnify catalogues')
    parser.add_argument('--output-accessions', 
                        '-o',
                        required=True,
                        type=Path,
                        help='File to write prepared list of non-redundant genome accessions')

    args = parser.parse_args()
    main(args.processed_acc, args.output_accessions, args.catalogue_metadata, args.gut_mapping)
