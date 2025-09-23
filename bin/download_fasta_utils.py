import csv
import logging
import shutil
from ftplib import FTP, error_perm, error_proto, error_reply, error_temp
from pathlib import Path

import boto3
import requests
from botocore import UNSIGNED
from botocore.config import Config
from botocore.exceptions import BotoCoreError, ClientError
from retry import retry


def get_fasta_url(accession: str, analysis_ftp_field: str = "generated_ftp") -> str:
    """
    Get the FTP URL for the fasta file of a given genome accession.
    For NCBI GCA accessions, return the NCBI datasets download link.
    For ENA accessions (ERZ or WGS), query the ENA API to get the FTP link.
    :param accession: Genome accession (ENA or NCBI)
    :param analysis_ftp_field: Name of the field containing FTP link (only for ERZ accessions). Can be one of:
        - generated_ftp
        - submitted_ftp
    :return: FTP URL of the fasta file
    """
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
        if file_url != "":
            return file_url
    raise ValueError(f"Empty string URL of the fasta file for accession {accession}")


@retry((BotoCoreError, ClientError), tries=5, delay=15, backoff=1.5)
def download_from_ENA_FIRE(accession: str, analysis_ftp_field: str, outpath: Path) -> Path:
    """
    Download the fasta file for a given ENA genome accession using the ENA FIRE S3.
    Only for ENA analysis accessions (ERZ).
    :param accession: ENA analysis accession
    :param analysis_ftp_field: Name of the field containing FTP link. Can be one of:
        - generated_ftp
        - submitted_ftp
    :param outpath: Path to the folder where the fasta file will be saved
    :return: Path to the downloaded fasta file
    """
    url = get_fasta_url(accession, analysis_ftp_field=analysis_ftp_field)
    logging.debug(f"Download {accession} from ENA FIRE using URL {url}")

    fire_endpoint = "http://hl.fire.sdo.ebi.ac.uk"
    fire_ena_bucket = "era-public"
    fire_path = url.replace("ftp.sra.ebi.ac.uk/vol1/", "")
    s3 = boto3.client("s3", endpoint_url=fire_endpoint, config=Config(signature_version=UNSIGNED))
    s3.download_file(fire_ena_bucket, fire_path, outpath)
    return check_if_empty_gz(outpath)


@retry((ValueError), tries=7, delay=15, backoff=2)
def download_from_ENA_API(accession: str, outpath: Path) -> Path:
    """
    Download the fasta file for a given genome accession using the ENA API.
    :param accession: accession
    :param outpath: Path to the folder where the fasta file will be saved
    :return: Path to the downloaded fasta file
    """
    api_endpoint = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}"
    logging.debug(f"Download {accession} from ENA API using URL {api_endpoint}")
    query = {"download": "true", "gzip": "true"}
    response = run_request(query, api_endpoint)

    with open(outpath, "wb") as out:
        out.write(response.content)
    return check_if_empty_gz(outpath)


@retry((OSError, error_temp, error_reply, error_proto, error_perm), tries=8, delay=10, backoff=3)
def download_from_ENA_FTP(accession: str, outpath: Path) -> Path:
    """
    Download the fasta file for a given ENA genome accession using the ENA FTP.
    :param accession: ENA accession
    :param outpath: Path to the folder where the fasta file will be saved
    :return: Path to the downloaded fasta file
    """
    url = get_fasta_url(accession)
    logging.debug(f"Download {accession} from ENA FTP using URL {url}")

    ftp_server = "ftp.ebi.ac.uk"
    ftp_path = url.replace(ftp_server, "")

    with FTP(ftp_server) as ftp:
        ftp.login()
        with open(outpath, "wb") as file:
            ftp.retrbinary(f"RETR {ftp_path}", file.write)
    return check_if_empty_gz(outpath)


def download_from_NCBI_datasets(accession: str, download_folder: Path) -> Path:
    """
    Download the assembly fasta file for a given NCBI genome accession using NCBI datasets.
    This function is a workaround, only to be used if other methods fail.
    :param accession: NCBI genome accession (GCA)
    :param download_folder: Path to the folder where the fasta file will be saved
    :return: Path to the downloaded fasta file
    """
    outpath = download_folder / f"{accession}.fa"
    accession_version = accession if "." in accession else accession + ".1"
    api_endpoint = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
        f"{accession_version}/download"
    )
    query = {"include_annotation_type": "GENOME_FASTA"}
    response = run_request(query, api_endpoint)

    content = response.read()
    tmp_archive_path = download_folder / "ncbi_tmp.zip"
    tmp_extract_path = download_folder / "ncbi_tmp"

    tmp_archive_path.write_bytes(content)

    shutil.unpack_archive(str(tmp_archive_path), str(tmp_extract_path))

    subdir_path = tmp_extract_path / f"ncbi_dataset/data/{accession_version}"
    source_files = list(subdir_path.glob("*_genomic.fna"))

    if not source_files:
        raise FileNotFoundError(f"No assembly fasta found in {subdir_path}")

    source_path = source_files[0]
    shutil.move(str(source_path), str(outpath))

    tmp_archive_path.unlink()
    shutil.rmtree(tmp_extract_path)

    if outpath.exists() and outpath.stat().st_size != 0:
        logging.debug(f"Successful. File saved to {outpath}")
        return outpath
    logging.debug(f"Downloaded file {outpath} has zero size. Removing the file.")
    outpath.unlink()
    raise ValueError(f"Downloaded file {outpath} has zero size")


def check_if_empty_gz(file_path: Path) -> Path:
    """
    Check if the downloaded gzipped fasta file is empty (less than 20 bytes size).
    If the file is empty, delete it and raise a ValueError.
    """
    if file_path.exists() and file_path.stat().st_size > 20:
        logging.debug(f"Successful. File saved to {file_path}")
        return file_path
    logging.debug(f"Downloaded file {file_path} has zero size. Removing the file.")
    file_path.unlink(missing_ok=True)
    raise ValueError(f"Downloaded file {file_path} has zero size.")


@retry(tries=5, delay=15, backoff=1.5)
def run_request(query, api_endpoint):
    """Run a GET request to the ENA search API with given query parameters."""
    response = requests.get(api_endpoint, params=query)
    response.raise_for_status()
    return response
