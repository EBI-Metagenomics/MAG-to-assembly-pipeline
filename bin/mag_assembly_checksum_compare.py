#!/usr/bin/env python3
# coding=utf-8

import argparse
import csv
import gzip 
import hashlib
import os
import requests
import shutil
import urllib.request as request
import urllib.parse

from Bio import SeqIO
import boto3
from botocore import UNSIGNED
from botocore.config import Config
from retry import retry
from tqdm import tqdm


def main(input_file, output_file, download_folder, minchecksum_match):
    completed_accessions = load_completed_accessions(output_file)
    with open(input_file, 'r') as inpt, open(output_file, 'a') as out:
        for line in tqdm(inpt.readlines()):
            line = line.strip().split('\t')
            mag_accession = line[0]
            mag_sample = line[1]
            derived_from_samples = line[2]
            assembly_accessions = line[3].split(",")
            if mag_accession in completed_accessions:
                continue

            if mag_accession.startswith("GCA"):
                # TODO add error exeptions to download functions
                if not "." in mag_accession:
                    mag_accession_version = mag_accession + ".1"
                else:
                    mag_accession_version = mag_accession
                mag_url = get_fasta_url(mag_accession_version)
                mag_file = download_fasta_from_ncbi(mag_url, download_folder, mag_accession_version) 
            else:
                mag_url = get_fasta_url(mag_accession)
                mag_file = download_fasta_from_ena(mag_url, download_folder, mag_accession)
            if mag_file:
                mag_hashes = compute_hashes(mag_file, write_cache=False)
                verified_assemblies = []
                for assembly_accession in assembly_accessions:
                    assembly_url = get_fasta_url(assembly_accession)
                    if not any(assembly_url.endswith(file_type) for file_type in [
                                                                                ".list.gz",
                                                                                ".vcf.gz",
                                                                                ".vcf",
                                                                                ".bam.gz",
                                                                                ".bam",
                                                                                ".bam.bai",
                                                                                ".bam.bai.gz",
                                                                                ".tsv",
                                                                                ".fastq.gz",
                                                                                ".fastq",
                                                                            ]):
                        assembly_file = download_fasta_from_ena(assembly_url, download_folder, assembly_accession)
                        if assembly_file:
                            assembly_hashes = compute_hashes(assembly_file, write_cache=True)
                            if is_assembly_correct(mag_hashes, assembly_hashes, minchecksum_match): 
                                verified_assemblies.append(assembly_accession)
                if verified_assemblies:
                    checksum_result = 'True - refined' if len(assembly_accessions) > 1 else 'True'
                    assembly_accessions = verified_assemblies
                else:
                    checksum_result = 'False'
                out.write('{}\t{}\t{}\t{}\t{}\n'.format(mag_accession, mag_sample, derived_from_samples, ",".join(assembly_accessions), checksum_result))


def load_completed_accessions(cache_path):
    try:
        with open(cache_path, "r") as f:
            return set(line.strip().split("\t")[0] for line in f.readlines())
    except FileNotFoundError:
        return set()


def is_fasta_file(file_path):
    try:
        with open(file_path, 'r') as handle:
            first_char = handle.read(1)
            return first_char == '>'
    except UnicodeDecodeError:
        os.remove(file_path)
        return False


@retry(tries=5, delay=15, backoff=1.5) 
def download_from_ENA_FIRE(url: str, download_folder: str, accession: str) -> str:
    outpath = os.path.join(download_folder, f'{accession}.fa.gz')
    cache_path = os.path.join(download_folder, f'{accession}.fa.hash')
    if (os.path.exists(outpath) and os.path.getsize(outpath) != 0) or \
        (os.path.exists(cache_path) and os.path.getsize(cache_path) != 0):
        return outpath
    
    if not os.path.exists(download_folder):
        os.makedirs(download_folder)

    fire_endpoint = "http://hl.fire.sdo.ebi.ac.uk"
    fire_ena_bucket = "era-public"
    fire_path = url.replace("ftp.sra.ebi.ac.uk/vol1/", "")

    # Create an S3 client with the Fire endpoint and unsigned requests
    s3 = boto3.client("s3", endpoint_url=fire_endpoint, config=Config(signature_version=UNSIGNED))
    
    # Download the file from the Fire storage
    s3.download_file(fire_ena_bucket, fire_path, outpath)

    if os.path.exists(outpath) and os.path.getsize(outpath) != 0:
        return outpath
    
    return None


@retry(tries=5, delay=15, backoff=1.5) 
def download_from_ENA_FTP(url, download_folder, accession):
    outpath = os.path.join(download_folder, f'{accession}.fa.gz')
    cache_path = os.path.join(download_folder, f'{accession}.fa.hash')
    if (os.path.exists(outpath) and os.path.getsize(outpath) != 0) or \
        (os.path.exists(cache_path) and os.path.getsize(cache_path) != 0):
        return outpath
    
    if not os.path.exists(download_folder):
        os.makedirs(download_folder)

    if not url.startswith('ftp://') or not url.startswith('https://'):
        url = 'https://' +  url

    response = request.urlopen(url)
    content = response.read()
    with open(outpath, 'wb') as out:
        out.write(content)

    if os.path.exists(outpath) and os.path.getsize(outpath) != 0:
            return outpath
    
    return None


@retry(tries=5, delay=15, backoff=1.5) 
def download_fasta_from_ncbi(url, download_folder, accession):
    outpath = os.path.join(download_folder, f'{accession}.fa')
    cache_path = os.path.join(download_folder, f'{accession}.fa.hash')
    if (os.path.exists(outpath) and os.path.getsize(outpath) != 0) or \
        (os.path.exists(cache_path) and os.path.getsize(cache_path) != 0):
        return outpath
    if not url.lower().startswith(('ftp', 'http')):
        print(url, 'is not an URL\n')
        return False
    if not os.path.exists(download_folder):
        os.makedirs(download_folder)
    response = request.urlopen(url)
    content = response.read()
    tmp_file = "ncbi_tmp.zip"
    tmp_path = os.path.join(download_folder, tmp_file)
    with open(tmp_path, 'wb') as out:
        out.write(content)
    shutil.unpack_archive(tmp_path, tmp_path.replace(".zip", ""))
    subdir_path = os.path.join(download_folder,"ncbi_tmp/ncbi_dataset/data/{}/".format(accession))
    source_file = [file for file in os.listdir(subdir_path) if file.endswith("_genomic.fna")]
    source_path = os.path.join(subdir_path, source_file[0]) # assembly_file is a list with one element
    shutil.move(source_path, outpath)
    # if not os.path.exists(outpath) or os.path.getsize(outpath) == 0:
    #     return None
    os.remove(tmp_path)
    shutil.rmtree(tmp_path.replace(".zip", ""))
    return outpath


def get_fasta_url(accession, analysis_ftp_field="generated_ftp"):
    if accession.startswith("GCA"):
        file_url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{}/download?include_annotation_type=GENOME_FASTA".format(accession)
        return file_url

    api_endpoint = 'https://www.ebi.ac.uk/ena/portal/api/search'
    accession_type = 'analysis' if accession.startswith("ERZ") else 'wgs_set'
    query = {
        "wgs_set": {
            "result": "wgs_set",
            "query": 'wgs_set="{}"'.format(accession),
            "fields": "set_fasta_ftp",
            "format": "tsv",
        },
        "analysis": {
            "result": "analysis",
            "query": 'analysis_accession="{}"'.format(accession),
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
    return None # no information about this accession in ENA


@retry(tries=5, delay=15, backoff=1.5)
def run_request(query, api_endpoint):
    request = requests.get(api_endpoint, params=urllib.parse.urlencode(query))
    request.raise_for_status()
    return request


def compute_hashes(file_path, write_cache=True, delete_fasta=True, separate_cache_dir=None):
    hashes = set()
    cache_path = file_path.replace(".fa.gz", ".fa") + ".hash"
    if separate_cache_dir:
        cache_basename = os.path.basename(cache_path)
        cache_path = os.path.join(separate_cache_dir, cache_basename)
    
    if os.path.exists(cache_path):
        with open(cache_path, "r") as handle:
            for line in handle:
                hashes.add(line.strip())
        return hashes

    if file_path.endswith(".gz"):
        with gzip.open(file_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                hash_object = hashlib.md5(str(record.seq.upper()).encode())
                hashes.add(hash_object.hexdigest())
    else:
        for record in SeqIO.parse(file_path, "fasta"):
            hash_object = hashlib.md5(str(record.seq.upper()).encode())
            hashes.add(hash_object.hexdigest())

    if write_cache:
        with open(cache_path, "w") as handle:
            for hash in hashes:
                handle.write(hash + "\n")

    if delete_fasta:
        os.remove(file_path)

    return hashes


def is_assembly_correct(mag_hash, assembly_hash, minchecksum_match):
    if minchecksum_match:
        # TODO what if MAG has less contigs than intersection size??
        intersection_size = len(mag_hash.intersection(assembly_hash))
        return intersection_size >= minchecksum_match
    return mag_hash.issubset(assembly_hash)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check correctness of assembly assignment for MAGs")
    parser.add_argument("input", help="Input TSV file with four columns: MAG accession, MAG sample accession, 'derived from' sample accession(s) and assembly accession(s) (both comma-separated)")
    parser.add_argument("output", help="Output TSV file with 5 columns: MAG accession, MAG sample accession, 'derived from' sample accession(s) and (refined) assembly accession and boolean result of checksum check")
    parser.add_argument("--download-folder", help="Folder to store downloaded files")
    parser.add_argument("--minchecksum-match", help="Minimun number of checksum matches between a MAG and an assembly to be correctly assigned. If stated, partial matching will be allowed", default=None, type=int)

    args = parser.parse_args()

    main(args.input, args.output, args.download_folder, args.minchecksum_match)
