import csv
import logging
import os
import shutil
import urllib.parse
from ftplib import FTP

import boto3
import requests
from botocore import UNSIGNED
from botocore.config import Config
from retry import retry


def get_fasta_url(accession, analysis_ftp_field="generated_ftp"):
    if accession.startswith("GCA"):
        file_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession}/download?include_annotation_type=GENOME_FASTA"
        return file_url

    api_endpoint = "https://www.ebi.ac.uk/ena/portal/api/search"
    accession_type = "analysis" if accession.startswith("ERZ") else "wgs_set"
    query = {
        "wgs_set": {
            "result": "wgs_set",
            "query": f'wgs_set="{accession}"',
            "fields": "set_fasta_ftp",
            "format": "tsv",
        },
        "analysis": {
            "result": "analysis",
            "query": f'analysis_accession="{accession}"',
            "fields": analysis_ftp_field,
            "format": "tsv",
        },
    }

    response = run_request(query[accession_type], api_endpoint)
    lines = response.text.splitlines()
    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        field_name = query[accession_type]["fields"]
        file_url = row[field_name].split(";")[0]  # Split to take the first FTP link if multiple
        return file_url
    return None  # no information about this accession in ENA


# TODO list errors explicitly, raise ValueError instead of returning None
@retry(tries=5, delay=15, backoff=1.5)
def download_from_ENA_FIRE(accession: str, analysis_ftp_field: str, outpath: str):
    url = get_fasta_url(accession, analysis_ftp_field=analysis_ftp_field)
    if not url:
        logging.debug(f"{accession} URL is empty for accession, ftp field: {analysis_ftp_field}")
        return None
        # raise ValueError(f"URL is empty, ftp field: {analysis_ftp_field}")
    logging.debug(f"Download {accession} from ENA FIRE using URL {url}")

    fire_endpoint = "http://hl.fire.sdo.ebi.ac.uk"
    fire_ena_bucket = "era-public"
    fire_path = url.replace("ftp.sra.ebi.ac.uk/vol1/", "")
    s3 = boto3.client("s3", endpoint_url=fire_endpoint, config=Config(signature_version=UNSIGNED))
    s3.download_file(fire_ena_bucket, fire_path, outpath)
    # 20 bytes is a size of an empty fa.gz
    if os.path.exists(outpath) and os.path.getsize(outpath) > 20:
        logging.debug(f"Successful. File saved to {outpath}")
        return outpath
    logging.debug(f"Downloaded file {outpath} has zero size. Removing the file.")
    os.remove(outpath)
    return None
    # raise ValueError(f"Downloaded file {outpath} has zero size")


@retry(tries=7, delay=15, backoff=2)
def download_from_ENA_API(accession: str, outpath: str) -> str:
    api_endpoint = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}"
    logging.debug(f"Download {accession} from ENA API using URL {api_endpoint}")
    query = {"download": "true", "gzip": "true"}
    response = requests.get(api_endpoint, params=urllib.parse.urlencode(query))
    response.raise_for_status()

    with open(outpath, "wb") as out:
        out.write(response.content)
    # 20 bytes is a size of an empty fa.gz
    if os.path.exists(outpath) and os.path.getsize(outpath) > 20:
        logging.debug(f"Successful. File saved to {outpath}")
        return outpath
    logging.debug(f"Downloaded file {outpath} has zero size. Removing the file.")
    os.remove(outpath)
    raise ValueError(f"Downloaded file {outpath} has zero size")


# TODO list errors explicitly, raise ValueError instead of returning None
@retry(tries=8, delay=10, backoff=3)
def download_from_ENA_FTP(accession, outpath):
    url = get_fasta_url(accession)
    if not url:
        logging.debug(f"{accession} URL is empty for accession")
        return None
        # raise ValueError(f"URL is empty")
    logging.debug(f"Download {accession} from ENA FTP using URL {url}")

    ftp_server = "ftp.ebi.ac.uk"
    ftp_path = url.replace(ftp_server, "")

    with FTP(ftp_server) as ftp:
        ftp.login()
        with open(outpath, "wb") as file:
            ftp.retrbinary(f"RETR {ftp_path}", file.write)
    # 20 bytes is a size of an empty fa.gz
    if os.path.exists(outpath) and os.path.getsize(outpath) > 20:
        logging.debug(f"Successful. File saved to {outpath}")
        return outpath
    logging.debug(f"Downloaded file {outpath} has zero size. Removing the file.")
    os.remove(outpath)
    return None
    # raise ValueError(f"Downloaded file {outpath} has zero size")


def download_from_NCBI_datasets(accession, download_folder):
    outpath = os.path.join(download_folder, f"{accession}.fa")
    accession_version = accession if "." in accession else accession + ".1"
    api_endpoint = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession_version}/download"
    query = {
        "include_annotation_type": "GENOME_FASTA",
    }
    response = run_request(query, api_endpoint)

    content = response.read()
    tmp_archive = "ncbi_tmp.zip"
    tmp_archive_path = os.path.join(download_folder, tmp_archive)
    tmp_path = tmp_archive_path.replace(".zip", "")
    with open(tmp_archive_path, "wb") as out:
        out.write(content)
    shutil.unpack_archive(tmp_archive_path, tmp_path)
    subdir_path = os.path.join(download_folder, f"ncbi_tmp/ncbi_dataset/data/{accession_version}/")
    source_file = [file for file in os.listdir(subdir_path) if file.endswith("_genomic.fna")]
    source_path = os.path.join(subdir_path, source_file[0])  # assembly_file is a list with one element
    shutil.move(source_path, outpath)
    os.remove(tmp_archive_path)
    shutil.rmtree(tmp_path)
    if os.path.exists(outpath) and os.path.getsize(outpath) != 0:
        logging.debug(f"Successful. File saved to {outpath}")
        return outpath
    logging.debug(f"Downloaded file {outpath} has zero size. Removing the file.")
    os.remove(outpath)
    raise ValueError(f"Downloaded file {outpath} has zero size")


@retry(tries=5, delay=15, backoff=1.5)
def run_request(query, api_endpoint):
    request = requests.get(api_endpoint, params=urllib.parse.urlencode(query))
    request.raise_for_status()
    return request
