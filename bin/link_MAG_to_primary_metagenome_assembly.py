#!/usr/bin/env python3
# coding=utf-8

import argparse
import json
import logging
import re
import shutil
import urllib.parse
from urllib.error import HTTPError

import requests
import xmltodict
from tqdm import tqdm
from retry import retry
from mag_assembly_checksum_compare import download_fasta_from_ena,download_fasta_from_ncbi,compute_hashes,get_fasta_url


# TODO maybe cache 'primary' or 'non-primary' type and run accessions for assemblies to avoid unnecessary API requests?


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        # logging.FileHandler(filename='script.log'),
        logging.StreamHandler()
    ]
)


def main(infile, outfile_confirmed, outfile_putative, outfile_fails, download_folder, cleanup, minchecksum_match):

    # load MAG accessions from existing out files to skip them
    completed_accessions = load_completed_accessions(outfile_confirmed, outfile_putative, outfile_fails, column_index=0)
    with open(infile, "r") as file_in, open(outfile_confirmed, "a") as out_confirmed, open(outfile_putative, "a") as out_putative, open(outfile_fails, 'a') as out_fails:
        for acc in tqdm(file_in.readlines()):     #  acc is a MAG accession 
            acc = acc.strip()
            if acc[:3] not in ["ERZ", "GCA"]:
                acc = acc.rstrip("0")
            if acc in completed_accessions:
                continue 

            # in this code block: obtain a corresponding sample id for the MAG id
            try:
                if acc.startswith("ERZ"):
                    mag_ena_data = load_data(acc, type="xml")
                    bin_sample = mag_ena_data['ANALYSIS_SET']['ANALYSIS']['SAMPLE_REF']['IDENTIFIERS']['PRIMARY_ID']
                else:
                    if acc.startswith("GCA"):   # if id is GenBank style it is required to convert it to ENA wgsSet id before requesting ENA API
                        mag_ena_data = load_data(genbank_to_ena_wgsset_accession(acc), type="summary")
                    else:     # remaining accessions are of wgsSet type
                        mag_ena_data = load_data(acc, type="summary")
                    bin_sample = mag_ena_data["summaries"][0]["sample"]
            except:
                print(acc, "unable to find sample accession", sep="\t", file=out_fails)
                continue

            # in this code block: for the determined sample id request upper level 'derived from' sample id(s) as well as run id(s) of the MAG
            try:
                sample_ena_data = load_data(bin_sample, type="xml")
                sample_attributes = sample_ena_data["SAMPLE_SET"]["SAMPLE"]["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]
            except:
                logging.info(f"Unable to find related samples for {acc}")
                print(acc, f"unable to load xml for {bin_sample}", sep="\t", file=out_fails)
                continue
            # attempting to retrieve either run id(s) or related sample id(s) from 'derived from' field
            derived_from_samples, derived_from_runs = extract_derived_from_info(sample_attributes)
            original_derived_from = derived_from_samples
            # this condition is required to find MAG samples that does not contain any needed information
            # or contain it in unexpected fields
            if not derived_from_samples and not derived_from_runs:
                print(acc, f"'derived from' field does not exist in XML or does not contain neither run nor sample accessions, MAG sample {bin_sample}", sep="\t", file=out_fails)
                continue   
            # if 'derived from' field does not contain any run id(s), look for them in the "description" field
            if not derived_from_runs:
                try:
                    description = sample_ena_data['SAMPLE_SET']['SAMPLE']['DESCRIPTION']
                    derived_from_runs = get_run_ids_from_description(description)
                    assert derived_from_runs
                except:
                    pass
                    # logging.info(f"Warning! Unable to identify runs which sample is derived from for MAG {acc}, {bin_sample}")
            # if "derived from" field does not contain any related sample id(s), look for them in the run(s) info using ENA API
            if not derived_from_samples:
                original_derived_from = derived_from_runs
                try:
                    derived_from_samples = get_samples_from_runs(derived_from_runs)
                    assert derived_from_samples
                except:
                    print(acc, f"unable to find 'derived from' sample from run accession for {bin_sample}", sep="\t", file=out_fails)
                    continue
            # in this code block: find all primary assemblies as well as their run ids linked to each of the related sample id(s)
            primary_assemblies, assembly2run = get_primary_metagenome_assembly_info(derived_from_samples)
            if not primary_assemblies:    # cases when primary assembly were not uploaded to ENA
                print(acc, f"there are no assemblies for sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}", sep="\t", file=out_fails)
                continue

            # in this code block: in cases when more than 1 assembly were found, try to descrease the number of assemblies
            # discarding samples which run id(s) does not match to MAG sample run id(s)
            if len(primary_assemblies) > 1 and derived_from_runs:
                try:
                    primary_assemblies = decrease_number_of_assemblies(primary_assemblies, derived_from_runs, assembly2run)
                except Exception as error:  # TODO improve this error handling
                    logging.info(f"Warning! Unable to decrease number of assembly accs for mag {acc}, sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}")
                    logging.info(str(error))

            # in this code block: verify retrieved assemblies using checksum comparason
            try:
                if acc.startswith("GCA"):
                    if not "." in acc:
                        acc_version = acc + '.1'
                    else:
                        acc_version = acc
                    mag_url = get_fasta_url(acc_version)
                    mag_file = download_fasta_from_ncbi(mag_url, download_folder, acc_version)
                else:
                    try:
                        # in ENA either generated_ftp or submitted_ftp (or both) fields may contain invalid links
                        mag_url = get_fasta_url(acc)
                        mag_file = download_fasta_from_ena(mag_url, download_folder, acc, unzip=True)
                    except HTTPError:
                        #retry downloading using submitted_ftp instead of generated_ftp field
                        mag_url = get_fasta_url(acc, analysis_ftp_field="submitted_ftp")
                        mag_file = download_fasta_from_ena(mag_url, download_folder, acc, unzip=True)
            except HTTPError as e:
                logging.info(f"HTTP Error while downloading MAG {acc}: {e.code} - {e.reason}")
                print(
                    acc, 
                    f"Failed to download MAG fasta file, sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}", 
                    sep="\t", file=out_fails
                )
                continue
            except Exception as e:
                logging.info(f"An error occurred during downloading of MAG {acc}:", e)
                print(
                    acc, 
                    f"Failed to download MAG fasta file, sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}",
                    sep="\t", file=out_fails
                )
                continue

            confirmed_assemblies = []
            putative_assemblies = {}
            mag_hashes = compute_hashes(mag_file, write_cash=False)
            for assembly in primary_assemblies:
                assembly_url = get_fasta_url(assembly)
                assembly_file = download_fasta_from_ena(assembly_url, download_folder, assembly, unzip=True)
                assembly_hashes = compute_hashes(assembly_file, write_cash=True)
                if mag_hashes.issubset(assembly_hashes): # TODO modify to avoid matching same MAG uploaded as ERZ as assembly as well as matching empty hashes
                    confirmed_assemblies.append(assembly)
                else:
                    intersection_size = len(mag_hashes.intersection(assembly_hashes))
                    if intersection_size >= minchecksum_match:
                        putative_assemblies[assembly] = intersection_size
                
            
            # Write a line to the output TSV file
            # Columns are      Genome_acc    Sample      Derived_from_sample     Derived_from_assembly
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

    if cleanup: 
        shutil.rmtree(download_folder)


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
            derived_from_samples.extend(re.findall("SAM[A-Z]+\d+|ERS\d+|SRS\d+|DRS\d+", attribute["VALUE"]))
            derived_from_runs.extend(re.findall("ERR\d+|SRR\d+|DRR\d+", attribute["VALUE"]))
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


def is_primary_metagenome(json_data):
    analysis_type = list(json_data["ANALYSIS_SET"]["ANALYSIS"]["ANALYSIS_TYPE"].keys())[0]
    if analysis_type == "SEQUENCE_ASSEMBLY":
        assembly_type = json_data["ANALYSIS_SET"]["ANALYSIS"]["ANALYSIS_TYPE"]["SEQUENCE_ASSEMBLY"]["TYPE"]
        if assembly_type == "primary metagenome":
            return True
    return False


def get_primary_metagenome_assembly_info(sample_accessions):
    primary_assemblies = set()
    assembly2run = {}
    api_endpoint = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    for sample_accession in sample_accessions:
        query = {
            'accession': '{}'.format(sample_accession),
            'result': 'analysis',
            'format': 'tsv'
        }
        response = run_request(query, api_endpoint)
        for line in response.text.splitlines():
            if not line.startswith("submitted"):
                line = line.strip().split("\t")
                submitted_ftp = line[0]
                analysis_accession = line[3]
                if not any(submitted_ftp.endswith(file_type) for file_type in [
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
                    assembly_data = load_data(analysis_accession, "xml")
                    if is_primary_metagenome(assembly_data):
                        runs = retrieve_assembly_runs_from_xml(assembly_data)
                        assembly2run[analysis_accession] = runs
                        primary_assemblies.add(analysis_accession)
    return primary_assemblies, assembly2run


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
                        help='Remove donloaded cash of checksums and download folder after execution')
    parser.add_argument("--minchecksum-match", default=0, type=int,
                        help="Minimun number of checksum matches between a MAG and an assembly to recognize as putative MAG-assembly pair. By default: 0")
    
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    main(args.infile, args.outfile_confirmed, args.outfile_putative, args.fails, args.download_folder, args.cleanup, args.minchecksum_match)
