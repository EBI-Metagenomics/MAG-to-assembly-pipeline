import logging 
from ftplib import FTP
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

work_dir = Path('/homes/sofia/nextflow/data')
catalogues_metadata_file = work_dir / 'all_catalog_metadata.tsv'
gut_mapping_file = work_dir / 'known_and_found_gut_genomes.tsv'
mag_layer_file = work_dir / "mag_layer.tsv"
bin_layer_file = work_dir / "bin_layer.tsv"
output_file = work_dir / "input_accessions.tsv"

previous = None


def main():
    logging.info("Starting preparation of the list of input accessions for genome-primary assembly linking script...")
    logging.info("Step 1/4:")
    logging.info("Download all genomes from MGnify catalogues...")
    # TODO what if nothing was downloaded??
    download_all_catalogues_metadata(work_dir, catalogues_metadata_file, previous)
    remove_gut_genomes(gut_mapping_file, catalogues_metadata_file)

    logging.info("Step 2/4:")
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

    logging.info("Step 3/4:")
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

    logging.info("Step 4/4:")
    logging.info("Merge all downloaded accessions removing redunduncy...")
    # TODO what if nothing was downloaded??   
    catalogues_df = pd.read_csv(catalogues_metadata_file, usecols=[2], sep='\t')
    mag_layer_df = pd.read_csv(mag_layer_file, usecols=[0,2], sep='\t')
    bin_layer_df = pd.read_csv(bin_layer_file, usecols=[0], sep='\t')

    ncbi_accessions = catalogues_df[catalogues_df['Genome_accession'].str.startswith('GCA')]['Genome_accession']
    mag_layer_df = mag_layer_df[~mag_layer_df['assembly_accession'].isin(ncbi_accessions)]

    combined_df = pd.concat([catalogues_df, mag_layer_df[['wgs_set']], bin_layer_df], ignore_index=True)
    combined_df = combined_df.drop_duplicates()
    combined_df.to_csv(output_file, sep='\t', index=False)

    logging.info(f"Finished successfully! List of genome accessions is saved to {output_file}")


def merge_without_redundancy():
    pass


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


# def download_layer(data, output_file):
#     url = "https://www.ebi.ac.uk/ena/portal/api/search"
#     try:
#         with requests.post(url, data=data, headers={"Content-Type": "application/x-www-form-urlencoded"}, stream=True) as response:
#             response.raise_for_status()  # Raise an error for bad status codes
#             with open(output_file, 'w') as file:
#                 for chunk in response.iter_content(chunk_size=1024):
#                     if chunk:
#                         file.write(chunk.decode('utf-8'))
#         logging.info(f"Data downloaded to {output_file}")
#     except requests.RequestException as e:
#         logging.error(f"Failed to download data: {e}")
    

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
    main()
