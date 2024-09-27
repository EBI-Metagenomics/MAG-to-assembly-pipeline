#!/usr/bin/env python3
# coding=utf-8

import argparse
import json
import logging
import os
import re
import shutil
import urllib.parse
from urllib.error import HTTPError

import requests
import xmltodict
from tqdm import tqdm
from retry import retry
from mag_assembly_checksum_compare import download_erz_from_ena,download_mag_from_ena,download_fasta_from_ncbi,compute_hashes,get_fasta_url


# TODO maybe cache run accessions for assemblies to avoid unnecessary API requests?

def setup_logging(debug=False):
    log_level = logging.DEBUG if debug else logging.INFO
    
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler()]
    )
    logging.getLogger('requests').setLevel(logging.WARNING)
    logging.getLogger('boto3').setLevel(logging.WARNING)
    logging.getLogger('botocore').setLevel(logging.WARNING)
    logging.getLogger('urllib').setLevel(logging.WARNING)
    logging.getLogger('urllib3').setLevel(logging.WARNING)
    logging.getLogger('s3transfer').setLevel(logging.WARNING)


def main(infile, outfile_confirmed, outfile_putative, outfile_fails, download_folder, cleanup, minchecksum_match):

    logging.debug("Loading MAG accessions from the input file as well as from existing output files to skip them...")
    completed_accessions = load_completed_accessions(outfile_confirmed, outfile_putative, outfile_fails, column_index=0)
    with open(infile, "r") as file_in, open(outfile_confirmed, "a") as out_confirmed, open(outfile_putative, "a") as out_putative, open(outfile_fails, 'a') as out_fails:
        for acc in tqdm(file_in.readlines()):     #  acc is a MAG accession 
            acc = acc.strip()
            logging.debug(f"Start processing of MAG/bin with accession {acc}")
            if acc[:3] not in ["ERZ", "GCA"]:
                acc = acc.rstrip("0")
            if acc in completed_accessions:
                logging.debug(f"The accession is in the list of completed accessions. Skipping")
                continue 

            logging.debug(f"Query ENA portal to get sample accession corresponding to the MAG/bin {acc}")
            try:
                if acc.startswith("ERZ"):
                    logging.debug(f"{acc} is an ENA analysis accession, retrieving metadata in XML from ENA portal")
                    mag_ena_data = load_data(acc, type="xml")
                    bin_sample = mag_ena_data['ANALYSIS_SET']['ANALYSIS']['SAMPLE_REF']['IDENTIFIERS']['PRIMARY_ID']
                elif acc.startswith("GCA"):
                    logging.debug(f"{acc} is a NCBI genome accession, retrieving summary for corresponding WGS set accession from ENA portal")
                    mag_ena_data = load_data(genbank_to_ena_wgsset_accession(acc), type="summary")
                else:
                    logging.debug(f"{acc} is a ENA WGS set accession, retrieving summary from ENA portal")
                    mag_ena_data = load_data(acc, type="summary")
                    bin_sample = mag_ena_data["summaries"][0]["sample"]
                logging.debug(f"Successful. MAG/bin {acc} sample accession is {bin_sample}")
            except:
                logging.debug(f"{acc} Unable to find sample accession. Skipping")
                print(acc, "unable to find sample accession", sep="\t", file=out_fails)
                continue

            # in this code block: for the determined sample id request upper level 'derived from' sample id(s) as well as run id(s) of the MAG
            logging.debug(f"Use ENA portal to find root sample accession and run accessions corresponding to the MAG/bin {acc}")
            try:
                logging.debug(f"Retrieving metadata in XML for accession {bin_sample} from ENA portal")
                sample_ena_data = load_data(bin_sample, type="xml")
                sample_attributes = sample_ena_data["SAMPLE_SET"]["SAMPLE"]["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]
            except:
                logging.info(f"Unable to find any related sample for {acc}. Skipping")
                print(acc, f"unable to load xml for {bin_sample}", sep="\t", file=out_fails)
                continue
            # attempting to retrieve either run id(s) or related sample id(s) from 'derived from' field
            logging.debug(f"Parsing sample attributes in XML metadata")
            derived_from_samples, derived_from_runs = extract_derived_from_info(sample_attributes)
            original_derived_from = derived_from_samples
            # this condition is required to find MAG samples that does not contain any needed information
            # or contain it in unexpected fields
            if not derived_from_samples and not derived_from_runs:
                print(acc, f"'derived from' field does not exist in XML or does not contain neither run nor sample accessions, MAG sample {bin_sample}", sep="\t", file=out_fails)
                logging.debug(f"'derived from' field does not exist in XML or does not contain neither run nor sample accessions, MAG sample {bin_sample}. Skipping")
                continue   
            logging.debug(f"According to the metadata bin sample was derived from {','.join(original_derived_from)}")
            # if 'derived from' field does not contain any run id(s), look for them in the "description" field
            if not derived_from_runs:
                logging.debug(f"'derived from' field does not contain any run id(s). Look for runs accessions in the metadata description")
                try:
                    description = sample_ena_data['SAMPLE_SET']['SAMPLE']['DESCRIPTION']
                    derived_from_runs = get_run_ids_from_description(description)
                    assert derived_from_runs
                    logging.debug(f"The following run accessions were found: {','.join(derived_from_runs)}")
                except:
                    logging.debug(f"Failed to identify run accessions for MAG/bin {acc}, {bin_sample}")
                    logging.debug("Comparson of run accessions for the MAG/bin {acc} and primary assemblies will be skipped")
            # if "derived from" field does not contain any related sample id(s), look for them in the run(s) info using ENA API
            if not derived_from_samples:
                logging.debug("'derived from' field does not contain sample accessions. Root sample will be identified through run accession")
                original_derived_from = derived_from_runs
                try:
                    derived_from_samples = get_samples_from_runs(derived_from_runs)
                    assert derived_from_samples
                    logging.debug(f"The following sample accessions were found: {','.join(derived_from_runs)}")
                except:
                    logging.debug(f"unable to find root sample from run accession for {bin_sample}. Skipping")
                    print(acc, f"unable to find 'derived from' sample from run accession for {bin_sample}", sep="\t", file=out_fails)
                    continue
            
            # in this code block: find all primary assemblies as well as their run ids linked to each of the related sample id(s)
            logging.debug(f"Find all primary metagenomic assemblies linked to the root sample {','.join(derived_from_runs)}")
            assembly2url = get_primary_assemblies_from_sample(derived_from_samples)
            primary_assemblies = list(assembly2url.keys())
            if not primary_assemblies:    # cases when primary assembly were not uploaded to ENA
                logging.debug(f"There are no assemblies for the given root sample")
                print(acc, f"there are no assemblies for sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}", sep="\t", file=out_fails)
                continue
            logging.debug(f"The following primary assemblies were found {','.join(primary_assemblies)}")

            # in this code block: in cases when more than 1 assembly were found, try to descrease the number of assemblies
            # discarding samples which run id(s) does not match to MAG sample run id(s)
            if len(primary_assemblies) > 1 and derived_from_runs:
                logging.debug(f"Attempt to decrease list of assemblies by filtering assemblies derived from the runs other than MAG runs.")
                assembly2run = get_run_accessions_for_assemblies(primary_assemblies)
                try:
                    primary_assemblies = decrease_number_of_assemblies(primary_assemblies, derived_from_runs, assembly2run)
                    logging.debug(f"Updated list of assemblies: {','.join(primary_assemblies)}")
                except Exception as error:  # TODO improve this error handling
                    logging.debug(f"Unable to decrease number of assemblies for MAG {acc}, sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}")
                    logging.debug(f"Due to {str(error)}")

            # in this code block: verify retrieved assemblies using checksum comparason
            logging.debug(f"Verify retrieved assemblies using comparason of contigs' hashes")
            try:
                if acc.startswith("GCA"):
                    if not "." in acc:
                        acc_version = acc + '.1'
                    else:
                        acc_version = acc
                    mag_url = get_fasta_url(acc_version)
                    logging.debug(f"Starting download of MAG/bin {acc} from NCBI using URL {mag_url}")
                    mag_file = download_fasta_from_ncbi(mag_url, download_folder, acc_version)
                    logging.debug(f"Successful. File saved to {mag_file}")
                else:
                    try:
                        # in ENA either generated_ftp or submitted_ftp (or both) fields may contain invalid links
                        mag_url = get_fasta_url(acc)
                        if acc.startswith("ERZ"):
                            logging.debug(f"Download MAG/bin {acc} from ENA FIRE using URL {mag_url}")
                            mag_file = download_erz_from_ena(mag_url, download_folder, acc)
                            logging.debug(f"Successful. File saved to {mag_file}")
                        else:
                            logging.debug(f"Download MAG/bin {acc} from ENA FTP using URL {mag_url}")
                            mag_file = download_mag_from_ena(mag_url, download_folder, acc)
                            logging.debug(f"Successful. File saved to {mag_file}")
                    except HTTPError:
                        logging.debug(f'Download from "generated_ftp" failed. Retry with "submitted_ftp"')
                        mag_url = get_fasta_url(acc, analysis_ftp_field="submitted_ftp")
                        if acc.startswith("ERZ"):
                            logging.debug(f"Download MAG/bin {acc} from ENA FIRE using URL {mag_url}")
                            mag_file = download_erz_from_ena(mag_url, download_folder, acc)
                            logging.debug(f"Successful. File saved to {mag_file}")
                        else:
                            logging.debug(f"Download MAG/bin {acc} from ENA FTP using URL {mag_url}")
                            mag_file = download_mag_from_ena(mag_url, download_folder, acc)
                            logging.debug(f"Successful. File saved to {mag_file}")
                    except botocore.exeptions.BotocoreError:
                        print("failed to download assembly", assembly, assembly_url)
                        continue

            except HTTPError as e:
                logging.info(f"HTTP Error while downloading MAG {acc}: {e.code} - {e.reason}. Skipping")
                print(
                    acc, 
                    f"Failed to download MAG fasta file, sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}", 
                    sep="\t", file=out_fails
                )
                continue
            except Exception as e:
                logging.info(f"An error occurred during downloading of MAG {acc}: {e}. Skipping")
                print(
                    acc, 
                    f"Failed to download MAG fasta file, sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}",
                    sep="\t", file=out_fails
                )
                continue

            confirmed_assemblies = []
            putative_assemblies = {}
            mag_hashes = compute_hashes(mag_file, write_cache=False)
            logging.debug(f"MAG/bin hashes were computed")
            for assembly in primary_assemblies:
                try:
                    assembly_url = assembly2url[assembly]
                    logging.debug(f"Download primary assembly {assembly} from ENA FTP using URL {assembly_url}")
                    assembly_file = download_erz_from_ena(assembly_url, download_folder, assembly)
                    logging.debug(f"Successful. File saved to {assembly_file}")
                except botocore.exceptions.ClientError:
                    print("failed to download assembly", assembly, assembly_url)
                    continue
                assembly_hashes = compute_hashes(assembly_file, write_cache=True)
                logging.debug(f"Assembly hashes were computed")
                if mag_hashes.issubset(assembly_hashes): # TODO modify to avoid matching empty file hashes
                    logging.debug(f"Assembly {assembly} is confirmed to be primary assembly for the MAG/bin {acc}")
                    confirmed_assemblies.append(assembly)
                else:
                    logging.debug(f"Assembly {assembly} is not a primary assembly for the MAG/bin {acc}")
                    intersection_size = len(mag_hashes.intersection(assembly_hashes))
                    if intersection_size >= minchecksum_match:
                        putative_assemblies[assembly] = intersection_size
                
            
            # Write a line to the output TSV file
            # Columns are      Genome_acc    Sample      Derived_from_sample     Derived_from_assembly
            logging.debug(f"Writing results to the output file")
            if confirmed_assemblies:
                print(
                    acc, 
                    bin_sample, 
                    ",".join(original_derived_from), 
                    ",".join(list(confirmed_assemblies)), 
                    sep="\t", file=out_confirmed
                )
            elif putative_assemblies:
                print(
                    acc,
                    bin_sample, 
                    ",".join(original_derived_from), 
                    ",".join(putative_assemblies),
                    sep="\t", file=out_putative
                )
            else:
                print(
                    acc, 
                    f"Found assemblies do not have sufficient matches, MAG sample {bin_sample}, 'derived from' sample {','.join(original_derived_from)}",
                    sep="\t", file=out_fails
                )
    
    if cleanup and os.path.exists(download_folder):
        shutil.rmtree(download_folder)
        logging.debug(f"Folder with downloaded files is deleted")


def load_completed_accessions(*files, column_index=0, filter_value="", separtor="\t"):
    completed_accessions = set()
    for file in files:
        try:
            with open(file, "r") as f:
                accessions = set(line.strip().split(separtor)[column_index] for line in f.readlines() if filter_value in line)
                completed_accessions.update(accessions)
        except FileNotFoundError:
            pass
    return completed_accessions


def genbank_to_ena_wgsset_accession(acc):
    try:
        summary_data = load_data(acc, type="summary")
        ena_accession = summary_data["summaries"][0]["wgsSet"]
        return ena_accession
    except:
       logging.info(f"Unable to convert NCBI accession {acc} to wgsSet accession")
       return None


def load_data(sample_id, type):
    url = 'https://www.ebi.ac.uk/ena/browser/api/{}/{}'.format(type, sample_id)
    try:
        request = run_browser_request(url)
    except:
        logging.info(f"Unable to request page content from URL for accession {sample_id}. Skipping.")
        return None
    if request.ok:
        try:
            if type == "xml":
                data_dict = xmltodict.parse(request.content)
                return json.loads(json.dumps(data_dict))
            elif type == "summary":
                return request.json()
        except:
            logging.info("Unable to load json from API for accession {}".format(sample_id))
            logging.info(request.text)
    else:
        logging.info('Could not retrieve xml for accession {}'.format(sample_id))
        logging.info(request.text)
        return None


def extract_derived_from_info(attributes):
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


def get_primary_assemblies_from_sample(sample_accessions):
    assembly2url = {}
    api_endpoint = "https://www.ebi.ac.uk/ena/portal/api/search"
    for sample_accession in sample_accessions:
        sample_type = "sample_accession" if sample_accession.startswith("SAM") else "secondary_sample_accession"
        query = {
            'result': 'analysis',
            'query': f'analysis_type=sequence_assembly AND assembly_type="primary metagenome" AND {sample_type}="{sample_accession}"',
            'format': 'tsv',
            'fields': 'generated_ftp,analysis_accession'
        }
        response = run_request(query, api_endpoint)
        for line in response.text.splitlines():
            if not line.startswith("generated_ftp"):
                assembly_url, assembly_accession = line.strip().split("\t")
                assembly_url = assembly_url.split(";")[0] 
                assembly2url[assembly_accession] = assembly_url
    return assembly2url


def get_run_accessions_for_assemblies(assembly_accessions):
    assembly2run = {}
    for assembly_accession in assembly_accessions:
        assembly_data = load_data(assembly_accession, "xml")
        runs = retrieve_assembly_runs_from_xml(assembly_data)
        assembly2run[assembly_accession] = runs
    return assembly2run


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


def decrease_number_of_assemblies(primary_assemblies, derived_from_runs, assembly2run):
    decreased_primary_assemblies = set()
    # logging.info("Trying to decrease number of putative assemblies for the MAG using run(s) accessions matching...")
    for assembly in primary_assemblies:
        assembly_run_acc = assembly2run[assembly]
        if not assembly_run_acc:
            decreased_primary_assemblies.add(assembly) # if it's impossible to find run(s) of the assembly, save it to check with checksum later
            continue
        if set(assembly_run_acc) == set(derived_from_runs):
            decreased_primary_assemblies.add(assembly)
    if decreased_primary_assemblies:
        return decreased_primary_assemblies
    return primary_assemblies


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


@retry(tries=5, delay=10, backoff=1.5)
def run_browser_request(url):
    request = requests.get(url)
    request.raise_for_status()
    return request


@retry(tries=5, delay=10, backoff=1.5)
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
