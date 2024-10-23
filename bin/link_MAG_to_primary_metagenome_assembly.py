#!/usr/bin/env python3
# coding=utf-8

import argparse
import csv
import json
import gzip
from ftplib import FTP
import logging
import os
import re
import shutil
import urllib.parse

import boto3
from botocore import UNSIGNED
from botocore.config import Config
from botocore.exceptions import ClientError, ParamValidationError
import requests
import xmltodict
from tqdm import tqdm
from retry import retry
from mag_assembly_checksum_compare import compute_hashes,get_fasta_url

# TODO add docs for functions and Type Annotations
# TODO look for primary assemblies even if bin sample is bio sample? 

def setup_logging(debug=False, error_logfile="ena_related_errors.log"):
    log_level = logging.DEBUG if debug else logging.INFO
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )
    
    error_handler = logging.FileHandler(error_logfile)
    error_handler.setLevel(logging.ERROR)
    simple_error_formatter = logging.Formatter('%(message)s')
    error_handler.setFormatter(simple_error_formatter)
    logging.getLogger().addHandler(error_handler)

    # Reduce logging for noisy libraries
    for noisy_lib in ['requests', 'boto3', 'botocore', 'urllib', 'urllib3', 's3transfer']:
        logging.getLogger(noisy_lib).setLevel(logging.WARNING)


def main(infile, outfile_confirmed, outfile_putative, outfile_fails, download_folder, cleanup, minchecksum_match):
    with open(infile, "r") as file_in, open(outfile_confirmed, "w") as out_confirmed, open(outfile_putative, "w") as out_putative, open(outfile_fails, 'w') as out_fails:
        for acc in tqdm(file_in.readlines()):     #  acc is a MAG/bin accession 
            acc = acc.strip() if acc[:3] in ["ERZ", "GCA"] else acc.strip().rstrip("0")

            logging.debug(f"Start processing of MAG/bin with accession {acc}")

            logging.debug(f"Query ENA portal to get bin sample accession corresponding to the MAG/bin {acc}")
            bin_sample = find_bin_sample_in_ena(acc)
            if not bin_sample:
                logging.info(f"{acc} Unable to find sample accession. Skipping")
                print(acc, "unable to find sample accession", sep="\t", file=out_fails)
                continue
            logging.debug(f"Successful. MAG/bin {acc} sample accession is {bin_sample}")
        
            logging.debug(f"Use ENA portal to find root sample accession and run accessions corresponding to the MAG/bin {acc}")
            derived_from, derived_from_samples, derived_from_runs = find_root_sample_and_run_in_ena(bin_sample)
            if not derived_from:
                print(acc, f"unable to load XML or 'derived from' field does not exist in XML, MAG sample {bin_sample}", sep="\t", file=out_fails)
                continue
            if not derived_from_samples:
                print(acc, f"unable to find 'derived from' sample from run metadata for {bin_sample}", sep="\t", file=out_fails)
                continue
            if not derived_from_runs:
                logging.debug(f"No bin's runs. Comparson of run accessions for the MAG/bin {acc} and primary assemblies will be skipped")
            
            logging.debug(f"Find all primary metagenomic assemblies linked to the root sample {','.join(derived_from_samples)}")
            primary_assemblies_dict = get_primary_assemblies_from_sample(derived_from_samples)
            if not primary_assemblies_dict:    # cases when primary assembly was not uploaded to ENA
                logging.debug(f"There are no assemblies for the given root sample")
                print(acc, f"there are no assemblies for sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}", sep="\t", file=out_fails)
                continue
            logging.debug(f"Successful. The following primary assemblies were found {','.join(primary_assemblies_dict.keys())}")

            if len(primary_assemblies_dict) > 1 and derived_from_runs:
                logging.debug(f"Attempt to decrease list of assemblies by filtering assemblies derived from the runs other than MAG runs")
                try:
                    primary_assemblies_dict = decrease_number_of_assemblies(primary_assemblies_dict, derived_from_runs)
                except Exception as error:  # TODO improve this error handling
                    logging.debug(f"Unable to decrease number of assemblies for MAG {acc}, sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}")
                    logging.debug(f"Due to {str(error)}")
                if not primary_assemblies_dict:
                    logging.info(f"All found primary assemblies were discarded during run comparason. Skipping")
                    print(acc, f"there are no assemblies with similar runs for sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}", sep="\t", file=out_fails)
                    continue
                logging.debug(f"Updated list of assemblies: {','.join(primary_assemblies_dict.keys())}")

            logging.debug(f"Verify retrieved assemblies using comparason of contigs' hashes")
            mag_hashes = handle_fasta_processing(acc, download_folder)
            if not mag_hashes:
                logging.info(f"Failed to download MAG {acc} fasta file. Skipping")
                print(
                    acc, 
                    f"Failed to download MAG fasta file, sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}",
                    sep="\t", file=out_fails
                )
                continue
            # continue
            logging.debug(f"MAG/bin hashes were computed")
            logging.debug(f"Start comparing MAG hashes to every primary assembly")
            confirmed_assemblies, putative_assemblies = compare_bin_and_assembly_hashes(acc, mag_hashes, primary_assemblies_dict, download_folder, minchecksum_match)
            logging.debug(f"Comparason finished")
            
            # Write a line to the output TSV file
            # Columns are      Genome_acc    Sample      Derived_from_sample     Derived_from_assembly
            logging.debug(f"Writing results to the output file")
            if confirmed_assemblies:
                print(acc, bin_sample, ",".join(derived_from), ",".join(confirmed_assemblies), sep="\t", file=out_confirmed)
            elif putative_assemblies:
                print(acc, bin_sample, ",".join(derived_from), ",".join(putative_assemblies), sep="\t", file=out_putative)
            else:
                print(acc, f"No sufficient matches for {bin_sample}, derived from {','.join(derived_from)}", file=out_fails)
    
    if cleanup and os.path.exists(download_folder):
        shutil.rmtree(download_folder)
        logging.debug(f"Folder with downloaded files is deleted")


def find_bin_sample_in_ena(acc):
    try:
        if acc.startswith("ERZ"):
            logging.debug(f"{acc} is an ENA analysis accession, retrieving metadata in XML from ENA portal")
            mag_ena_data = load_data(acc, type="xml")
            return mag_ena_data['ANALYSIS_SET']['ANALYSIS']['SAMPLE_REF']['IDENTIFIERS']['PRIMARY_ID']
        elif acc.startswith("GCA"):
            logging.debug(f"{acc} is a NCBI genome accession, retrieving metadata in XML from ENA portal")
            mag_ena_data = load_data(acc, type="xml")
            return mag_ena_data['ASSEMBLY_SET']['ASSEMBLY']['SAMPLE_REF']['IDENTIFIERS']['PRIMARY_ID']
        else:
            logging.debug(f"{acc} is an ENA WGS set accession, retrieving summary from ENA portal")
            mag_ena_data = load_data(acc, type="summary")
            return mag_ena_data["summaries"][0]["sample"]
    except Exception as e:
        logging.debug(f"Failed to fetch sample accession for {acc} due to {e}")
        return None


def find_root_sample_and_run_in_ena(bin_sample):
    try:
        logging.debug(f"Retrieving metadata in XML for accession {bin_sample} from ENA portal")
        sample_ena_data = load_data(bin_sample, type="xml")
        sample_attributes = sample_ena_data["SAMPLE_SET"]["SAMPLE"]["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]
        logging.debug(f"Parsing sample attributes in XML metadata")
        derived_from_samples, derived_from_runs = parse_derived_from_attribute(sample_attributes)
        assert derived_from_runs or derived_from_samples, "No 'derived from' attribute"
    except AssertionError as e:
        logging.debug(f"Unable to parse sample XML attributes for {bin_sample} due to: {e}")
        return None, None, None
    except Exception as e:
        logging.info(f"Unable to get bin sample XML or parse its attributes for {bin_sample} due to: {e}")
        return None, None, None
    
    derived_from = derived_from_samples if derived_from_samples else derived_from_runs
    logging.debug(f"According to the metadata bin sample was derived from {','.join(derived_from)}")

    # if 'derived from' field does not contain any run id(s), look for them in the "description" field
    if not derived_from_runs:
        logging.debug(f"Look for runs accessions in the bin sample metadata <DESCRIPTION> field")
        try:
            description = sample_ena_data['SAMPLE_SET']['SAMPLE']['DESCRIPTION']
            derived_from_runs = get_run_ids_from_description(description)
            assert derived_from_runs
            logging.debug(f"The following run accessions were found: {','.join(derived_from_runs)}")
            return derived_from, derived_from_samples, derived_from_runs
        except:
            logging.debug(f"Failed to identify run accessions for bin sample {bin_sample}")
            return derived_from, derived_from_samples, None
    
    # if "derived from" field does not contain any related sample id(s), look for them in the run(s) metadata using ENA API
    if not derived_from_samples:
        logging.debug("Root sample will be identified through run accession(s)")
        try:
            derived_from_samples = get_samples_from_runs(derived_from_runs)
            assert derived_from_samples
            logging.debug(f"The following sample accessions were found: {','.join(derived_from_samples)}")
            return derived_from, derived_from_samples, derived_from_runs
        except:
            logging.debug(f"unable to find root sample from run accession for {bin_sample}")
            return derived_from, None, derived_from_runs
        
    return derived_from, derived_from_samples, derived_from_runs


def get_primary_assemblies_from_sample(sample_accessions):
    primary_assemblies_dict = {}
    api_endpoint = "https://www.ebi.ac.uk/ena/portal/api/search"
    for sample_accession in sample_accessions:
        sample_type = "sample_accession" if sample_accession.startswith("SAM") else "secondary_sample_accession"
        query = {
            'result': 'analysis',
            'query': f'analysis_type=sequence_assembly AND assembly_type="primary metagenome" AND {sample_type}="{sample_accession}"',
            'format': 'tsv',
            'fields': 'generated_ftp,run_accession,analysis_accession'
        }
        response = run_request(query, api_endpoint)
        lines = response.text.splitlines()
        reader = csv.DictReader(lines, delimiter="\t")
        
        for row in reader:
            assembly_url = row['generated_ftp'].split(";")[0]  # Split to take the first FTP link if multiple
            run_accession = row['run_accession']
            assembly_accession = row['analysis_accession']
            primary_assemblies_dict[assembly_accession] = (assembly_url, run_accession)
    
    return primary_assemblies_dict


def retrieve_assembly_runs_from_xml(assembly_data):
    try:
        run_ref_data = assembly_data["ANALYSIS_SET"]["ANALYSIS"]["RUN_REF"]
        if isinstance(run_ref_data, list):
            return [run["IDENTIFIERS"]["PRIMARY_ID"] for run in run_ref_data]
        else: 
            return [run_ref_data["IDENTIFIERS"]["PRIMARY_ID"]]
    except KeyError:
        try:
            analysis_description = assembly_data['ANALYSIS_SET']['ANALYSIS']['DESCRIPTION']
            return get_run_ids_from_description(analysis_description)
        except:
            return None


def decrease_number_of_assemblies(assembly2metadata, bin_runs):
    for assembly, (_, assembly_run_acc) in list(assembly2metadata.items()):
        if not assembly_run_acc: # if run(s) of the assembly not found, save it to check with checksum later
            continue
        if {assembly_run_acc} != set(bin_runs):
            del assembly2metadata[assembly]

    return assembly2metadata


def handle_fasta_processing(accession, download_folder):
    try:
        outpath = os.path.join(download_folder, f'{accession}.fa.gz')
        cache_path = os.path.join(download_folder, f'{accession}.fa.hash')
        if (os.path.exists(outpath) and os.path.getsize(outpath) != 0) or \
            (os.path.exists(cache_path) and os.path.getsize(cache_path) != 0):
            return compute_hashes(outpath, write_cache=False)
        
        if not os.path.exists(download_folder):
            os.makedirs(download_folder)
            logging.debug(f"Directory {download_folder} is created")

        if accession.startswith("ERZ"):
            # in ENA generated_ftp or submitted_ftp (or both) fields may contain invalid links
            try:
                fasta_file = download_from_ENA_FIRE(accession, "generated_ftp", outpath)
                if fasta_file is None:
                    raise ValueError("Empty URL or empty file in 'generated_ftp'")
                return compute_hashes(fasta_file, write_cache=True)
            except (gzip.BadGzipFile, ClientError, ParamValidationError, ValueError) as e:
                logging.error(f"{accession} Download from link in 'generated_ftp' failed due to: {e}")
                logging.debug(f'Retry with "submitted_ftp"')
                fasta_file = download_from_ENA_FIRE(accession, "submitted_ftp", outpath)
                if fasta_file is None:
                    raise ValueError("Empty URL or empty file in 'submitted_ftp'")
                return compute_hashes(fasta_file, write_cache=True)
        elif accession.startswith("GCA"):
            fasta_file = download_from_ENA_API(accession, outpath)
            return compute_hashes(fasta_file, write_cache=False)
        else:
            fasta_file = download_from_ENA_FTP(accession, outpath)
            if fasta_file is None:
                raise ValueError("Empty URL or empty file'")
            return compute_hashes(fasta_file, write_cache=False)

    except requests.HTTPError as e:
        logging.error(f"{accession} HTTP Error while downloading: {e.code} - {e.reason}")
        return None
    except Exception as e:
        logging.error(f"{accession} Failed to process fasta file due to: {e}")
        return None


@retry(tries=5, delay=15, backoff=1.5) 
def download_from_ENA_FIRE(accession: str, analysis_ftp_field: str, outpath: str) -> str:
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
    query = {
        'download': 'true',
        'gzip': 'true'
    }
    response = requests.get(api_endpoint, params=urllib.parse.urlencode(query))
    response.raise_for_status()
    
    with open(outpath, 'wb') as out:
        out.write(response.content)
    # 20 bytes is a size of an empty fa.gz
    if os.path.exists(outpath) and os.path.getsize(outpath) > 20:
        logging.debug(f"Successful. File saved to {outpath}")
        return outpath
    logging.debug(f"Downloaded file {outpath} has zero size. Removing the file.")
    os.remove(outpath)
    raise ValueError(f"Downloaded file {outpath} has zero size")


def download_from_NCBI_datasets(accession, download_folder):
    outpath = os.path.join(download_folder, f'{accession}.fa')
    accession_version = accession if "." in accession else accession + ".1"
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{accession_version}/download"
    query = {
        'include_annotation_type': 'GENOME_FASTA',
    }
    response = requests.get(url, params=urllib.parse.urlencode(query))
    response.raise_for_status()

    content = response.read()
    tmp_archive = "ncbi_tmp.zip"
    tmp_archive_path = os.path.join(download_folder, tmp_archive)
    tmp_path = tmp_archive_path.replace(".zip", "")
    with open(tmp_archive_path, 'wb') as out:
        out.write(content)
    shutil.unpack_archive(tmp_archive_path, tmp_path)
    subdir_path = os.path.join(download_folder, f"ncbi_tmp/ncbi_dataset/data/{accession_version}/")
    source_file = [file for file in os.listdir(subdir_path) if file.endswith("_genomic.fna")]
    source_path = os.path.join(subdir_path, source_file[0]) # assembly_file is a list with one element
    shutil.move(source_path, outpath)
    os.remove(tmp_archive_path)
    shutil.rmtree(tmp_path)
    if os.path.exists(outpath) and os.path.getsize(outpath) != 0:
        logging.debug(f"Successful. File saved to {outpath}")
        return outpath
    logging.debug(f"Downloaded file {outpath} has zero size. Removing the file.")
    os.remove(outpath)
    raise ValueError(f"Downloaded file {outpath} has zero size")


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
        with open(outpath, 'wb') as file:
            ftp.retrbinary(f"RETR {ftp_path}", file.write)
    # 20 bytes is a size of an empty fa.gz
    if os.path.exists(outpath) and os.path.getsize(outpath) > 20:
        logging.debug(f"Successful. File saved to {outpath}")
        return outpath
    logging.debug(f"Downloaded file {outpath} has zero size. Removing the file.")
    os.remove(outpath)
    return None
    # raise ValueError(f"Downloaded file {outpath} has zero size")


def compare_bin_and_assembly_hashes(acc, mag_hashes, assembly2metadata, download_folder, minchecksum_match):
    confirmed_assemblies = []
    putative_assemblies = {}
    for assembly in assembly2metadata:
        assembly_hashes = handle_fasta_processing(assembly, download_folder)
        if not assembly_hashes:
            logging.info(f"For the MAG {acc} failed to download primary assembly {assembly} fasta file.")
            continue
        logging.debug(f"Assembly hashes were computed")
        if mag_hashes.issubset(assembly_hashes): # TODO modify to avoid matching empty file hashes
            logging.debug(f"Assembly {assembly} is confirmed to be primary assembly for the MAG/bin {acc}")
            confirmed_assemblies.append(assembly)
        else:
            logging.debug(f"Assembly {assembly} is not a primary assembly for the MAG/bin {acc}")
            intersection_size = len(mag_hashes.intersection(assembly_hashes))
            if intersection_size >= minchecksum_match:
                putative_assemblies[assembly] = intersection_size
    return confirmed_assemblies, putative_assemblies


def genbank_to_ena_wgsset_accession(acc):
    try:
        summary_data = load_data(acc, type="summary")
        ena_accession = summary_data["summaries"][0]["wgsSet"]
        return ena_accession
    except:
       logging.info(f"Unable to convert NCBI accession {acc} to wgsSet accession")
       return None


def load_data(accession, type):
    url = f'https://www.ebi.ac.uk/ena/browser/api/{type}/{accession}'
    try:
        request = run_browser_request(url)
        if type == "xml":
                data_dict = xmltodict.parse(request.content)
                return json.loads(json.dumps(data_dict))
        elif type == "summary":
                return request.json()
    except Exception as e:
        logging.error(f"{accession} Unable to request page content from URL {url} due to: {e}")
        return None


def parse_derived_from_attribute(attributes):
    derived_from_samples = []
    derived_from_runs = []

    for attribute in attributes:
        if all([x in attribute["TAG"] for x in ["derived", "from"]]):
            derived_from_samples.extend(re.findall(r"SAM[A-Z]+\d+|ERS\d+|SRS\d+|DRS\d+", attribute["VALUE"]))
            derived_from_runs.extend(re.findall(r"ERR\d+|SRR\d+|DRR\d+", attribute["VALUE"]))
            break

    return derived_from_samples, derived_from_runs


def get_samples_from_runs(runs):
    samples = []
    for run in runs:
        run_data = load_data(run, "xml")
        for link in run_data["RUN_SET"]["RUN"]["RUN_LINKS"]["RUN_LINK"]:
            if link["XREF_LINK"]["DB"] == "ENA-SAMPLE":
                sample = link["XREF_LINK"]["ID"]
                samples.append(sample)
                break
    return samples


def get_run_ids_from_description(description):
    def unfold_accession_range(start, end):
        start_num = int(start[3:])  # Extract the numeric part after the prefix (e.g., ERR)
        end_num = int(end[3:])
        return [f"{start[:3]}{num}" for num in range(start_num, end_num + 1)]

    matches = re.findall(r'\b(?:ERR|SRR|DRR)\d+(?:-(?:ERR|SRR|DRR)\d+)?|\b(?:ERR|SRR|DRR)\d+\b', description)
    unfolded_accessions = []
    for match in matches:
        if '-' in match:  # If the match is a range
            start, end = match.split('-')
            unfolded_accessions.extend(unfold_accession_range(start, end))
        else:
            unfolded_accessions.append(match)

    return unfolded_accessions


@retry(tries=5, delay=15, backoff=1.5)
def run_browser_request(url):
    request = requests.get(url)
    request.raise_for_status()
    return request


@retry(tries=5, delay=15, backoff=1.5)
def run_request(query, api_endpoint):
    request = requests.get(api_endpoint, params=urllib.parse.urlencode(query))
    request.raise_for_status()
    return request


def parse_args():
    parser = argparse.ArgumentParser(description='The script takes in a file with bin/MAG accessions'
                                                 '(ERZ, GCA, wgsSet) and outputs the accession of the '
                                                 'primary metagenome each came from where possible')
    parser.add_argument('-i', '--infile', required=True,
                        help='Path to the file containing a list of bin/MAG accessions, one per line.')
    parser.add_argument('-o', '--outfile-confirmed', required=True,
                        help='Name of the outfile to write MAG-assembly pairs that have full checksum matching. If file exists, lines will be added to the end of the file')
    parser.add_argument('-p', '--outfile-putative', required=True,
                        help='Name of the outfile to write MAG-assembly pairs that have partial checksum matching. If file exists, lines will be added to the end of the file')
    parser.add_argument('-f', '--fails', default='failed_accessions.tsv',
                        help='File to output failed MAG accessions.  If file exists, lines will be added to the end of the file. By default: failed_accessions.tsv')
    parser.add_argument("--download-folder",required=True, 
                        help='Folder to store downloaded files. By default: fasta_downloads', 
                        default='fasta_downloads')
    parser.add_argument("--cleanup", action="store_true", 
                        help='Remove downloaded cache of checksums and download folder after execution')
    parser.add_argument("--minchecksum-match", default=0, type=int,
                        help="Minimun number of checksum matches between a MAG and an assembly to recognize as putative MAG-assembly pair. By default: 0")
    parser.add_argument("--debug", action="store_true", 
                        help='Print out more information')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    setup_logging(args.debug)
    main(args.infile, args.outfile_confirmed, args.outfile_putative, args.fails, args.download_folder, args.cleanup, args.minchecksum_match)
