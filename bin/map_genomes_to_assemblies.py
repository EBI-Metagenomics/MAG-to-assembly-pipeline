#!/usr/bin/env python3
# coding=utf-8

import argparse
import csv
import json
import logging
import re
import urllib.parse

import requests
import xmltodict
from retry import retry
from tqdm import tqdm

# TODO add docs for functions and Type Annotations
# TODO look for primary assemblies even if bin sample is bio sample?


def setup_logging(debug=False, error_logfile="ena_related_errors.log"):
    log_level = logging.DEBUG if debug else logging.INFO

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler()],
    )

    error_handler = logging.FileHandler(error_logfile)
    error_handler.setLevel(logging.ERROR)
    simple_error_formatter = logging.Formatter("%(message)s")
    error_handler.setFormatter(simple_error_formatter)
    logging.getLogger().addHandler(error_handler)

    # Reduce logging for noisy libraries
    for noisy_lib in ["requests", "urllib", "urllib3"]:
        logging.getLogger(noisy_lib).setLevel(logging.WARNING)


def main(infile, outfile):
    err_lines = []
    out_lines = []
    with open(infile, "r") as file_in:
        for genome_accession in tqdm(file_in.readlines()):
            genome_accession = genome_accession.strip()
            if genome_accession[:3] not in ["ERZ", "GCA"]:
                genome_accession = genome_accession.rstrip("0")  # CAMPAA010000000 -> CAMPAA01

            logging.debug(f"Start processing of MAG/bin with accession {genome_accession}")

            logging.debug(
                f"Query ENA API to get bin sample accession corresponding to the MAG/bin {genome_accession}"
            )
            bin_sample = find_bin_sample_in_ena(genome_accession)
            if not bin_sample:
                logging.info(f"{genome_accession} Unable to find bin sample accession. Skipping")
                err_lines.append(f"{genome_accession}\tunable to find sample accession")
                continue
            logging.debug(
                f"Successful. MAG/bin {genome_accession} sample accession is {bin_sample}"
            )

            logging.debug(
                f"Use ENA API to find root sample accession and run accessions corresponding to the MAG/bin {genome_accession}"
            )
            derived_from, derived_from_samples, derived_from_runs = find_root_sample_and_run_in_ena(
                bin_sample
            )
            if not derived_from:
                err_lines.append(
                    f"{genome_accession}\tunable to load XML or 'derived from' field does not exist in XML, MAG sample {bin_sample}"
                )
                continue
            if not derived_from_samples:
                err_lines.append(
                    f"{genome_accession}\tunable to find 'derived from' sample from run metadata for {bin_sample}"
                )
                continue
            if not derived_from_runs:
                logging.debug(
                    f"No bin's runs. Comparson of run accessions for the MAG/bin {genome_accession} and primary assemblies will be skipped"
                )

            logging.debug(
                f"Find all primary metagenomic assemblies linked to the root sample {','.join(derived_from_samples)}"
            )
            primary_assemblies = get_primary_assemblies_from_sample(derived_from_samples)
            if not primary_assemblies:  # cases when primary assembly was not uploaded to ENA
                logging.debug("There are no assemblies for the given root sample")
                err_lines.append(
                    f"{genome_accession}\tthere are no assemblies for sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}"
                )
                continue
            logging.debug(
                f"Successful. The following primary assemblies were found {','.join(primary_assemblies)}"
            )

            if len(primary_assemblies) > 1 and derived_from_runs:
                logging.debug(
                    "Attempt to decrease list of assemblies by filtering assemblies derived from the runs other than MAG runs"
                )
                try:
                    primary_assemblies = decrease_number_of_assemblies(
                        primary_assemblies, derived_from_runs
                    )
                except Exception as error:  # TODO improve this error handling
                    logging.debug(
                        f"Unable to decrease number of assemblies for MAG {genome_accession}, sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}"
                    )
                    logging.debug(f"Due to {str(error)}")
                if not primary_assemblies:
                    logging.info(
                        "All found primary assemblies were discarded during run comparason. Skipping"
                    )
                    err_lines.append(
                        f"{genome_accession}\tthere are no assemblies with similar runs for sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}"
                    )
                    continue
                logging.debug(f"Updated list of assemblies: {','.join(primary_assemblies)}")

            if primary_assemblies:
                logging.debug("Write list of assemblies to the output file")
                out_lines.append(f"{genome_accession}\t{','.join(primary_assemblies)}")
            else:
                logging.debug(f"No primary assemblies found for {genome_accession}")
                err_lines.append(
                    f"{genome_accession}\tno primary assemblies found for sample id: {bin_sample}, derived samples: {','.join(derived_from_samples)}"
                )


def find_bin_sample_in_ena(genome_accession):
    try:
        if genome_accession.startswith("ERZ"):
            logging.debug(
                f"{genome_accession} is an ENA analysis accession, retrieving metadata in XML from ENA portal"
            )
            mag_ena_data = load_data(genome_accession, type="xml")
            return mag_ena_data["ANALYSIS_SET"]["ANALYSIS"]["SAMPLE_REF"]["IDENTIFIERS"][
                "PRIMARY_ID"
            ]
        elif genome_accession.startswith("GCA"):
            logging.debug(
                f"{genome_accession} is a NCBI genome accession, retrieving metadata in XML from ENA portal"
            )
            mag_ena_data = load_data(genome_accession, type="xml")
            return mag_ena_data["ASSEMBLY_SET"]["ASSEMBLY"]["SAMPLE_REF"]["IDENTIFIERS"][
                "PRIMARY_ID"
            ]
        else:
            logging.debug(
                f"{genome_accession} is an ENA WGS set accession, retrieving summary from ENA portal"
            )
            mag_ena_data = load_data(genome_accession, type="summary")
            return mag_ena_data["summaries"][0]["sample"]
    except Exception as e:
        logging.debug(f"Failed to fetch sample accession for {genome_accession} due to {e}")
        return None


def find_root_sample_and_run_in_ena(bin_sample):
    try:
        logging.debug(f"Retrieving metadata in XML for accession {bin_sample} from ENA portal")
        sample_ena_data = load_data(bin_sample, type="xml")
        sample_attributes = sample_ena_data["SAMPLE_SET"]["SAMPLE"]["SAMPLE_ATTRIBUTES"][
            "SAMPLE_ATTRIBUTE"
        ]
        logging.debug("Parsing sample attributes in XML metadata")
        derived_from_samples, derived_from_runs = parse_derived_from_attribute(sample_attributes)
        assert derived_from_runs or derived_from_samples, "No 'derived from' attribute"
    except AssertionError as e:
        logging.debug(f"Unable to parse sample XML attributes for {bin_sample} due to: {e}")
        return None, None, None
    except Exception as e:
        logging.info(
            f"Unable to get bin sample XML or parse its attributes for {bin_sample} due to: {e}"
        )
        return None, None, None

    derived_from = derived_from_samples if derived_from_samples else derived_from_runs
    logging.debug(
        f"genome_accessionording to the metadata bin sample was derived from {','.join(derived_from)}"
    )

    # if 'derived from' field does not contain any run id(s), look for them in the "description" field
    if not derived_from_runs:
        logging.debug("Look for runs accessions in the bin sample metadata <DESCRIPTION> field")
        try:
            description = sample_ena_data["SAMPLE_SET"]["SAMPLE"]["DESCRIPTION"]
            derived_from_runs = get_run_ids_from_description(description)
            assert derived_from_runs
            logging.debug(f"The following run accessions were found: {','.join(derived_from_runs)}")
            return derived_from, derived_from_samples, derived_from_runs
        except Exception:
            logging.debug(f"Failed to identify run accessions for bin sample {bin_sample}")
            return derived_from, derived_from_samples, None

    # if "derived from" field does not contain any related sample id(s), look for them in the run(s) metadata using ENA API
    if not derived_from_samples:
        logging.debug("Root sample will be identified through run accession(s)")
        try:
            derived_from_samples = get_samples_from_runs(derived_from_runs)
            assert derived_from_samples
            logging.debug(
                f"The following sample accessions were found: {','.join(derived_from_samples)}"
            )
            return derived_from, derived_from_samples, derived_from_runs
        except Exception:
            logging.debug(f"unable to find root sample from run accession for {bin_sample}")
            return derived_from, None, derived_from_runs

    return derived_from, derived_from_samples, derived_from_runs


def get_primary_assemblies_from_sample(sample_accessions):
    primary_assemblies = []
    api_endpoint = "https://www.ebi.ac.uk/ena/portal/api/search"
    for sample_accession in sample_accessions:
        sample_type = (
            "sample_accession"
            if sample_accession.startswith("SAM")
            else "secondary_sample_accession"
        )
        query = {
            "result": "analysis",
            "query": f'analysis_type=sequence_assembly AND assembly_type="primary metagenome" AND {sample_type}="{sample_accession}"',
            "format": "tsv",
            "fields": "generated_ftp,run_accession,analysis_accession",
        }
        response = run_request(query, api_endpoint)
        lines = response.text.splitlines()
        reader = csv.DictReader(lines, delimiter="\t")

        for row in reader:
            assembly_accession = row["analysis_accession"]
            primary_assemblies.append(assembly_accession)

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
            analysis_description = assembly_data["ANALYSIS_SET"]["ANALYSIS"]["DESCRIPTION"]
            return get_run_ids_from_description(analysis_description)
        except Exception:
            return None


def decrease_number_of_assemblies(assembly2metadata, bin_runs):
    for assembly, (_, assembly_run_genome_accession) in list(assembly2metadata.items()):
        if (
            not assembly_run_genome_accession
        ):  # if run(s) of the assembly not found, save it to check with checksum later
            continue
        if {assembly_run_genome_accession} != set(bin_runs):
            del assembly2metadata[assembly]

    return assembly2metadata


def genbank_to_ena_wgsset_accession(genome_accession):
    try:
        summary_data = load_data(genome_accession, type="summary")
        ena_accession = summary_data["summaries"][0]["wgsSet"]
        return ena_accession
    except Exception:
        logging.info(f"Unable to convert NCBI accession {genome_accession} to wgsSet accession")
        return None


def load_data(accession, type):
    url = f"https://www.ebi.ac.uk/ena/browser/api/{type}/{accession}"
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
            derived_from_samples.extend(
                re.findall(r"SAM[A-Z]+\d+|ERS\d+|SRS\d+|DRS\d+", attribute["VALUE"])
            )
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

    matches = re.findall(
        r"\b(?:ERR|SRR|DRR)\d+(?:-(?:ERR|SRR|DRR)\d+)?|\b(?:ERR|SRR|DRR)\d+\b", description
    )
    unfolded_accessions = []
    for match in matches:
        if "-" in match:  # If the match is a range
            start, end = match.split("-")
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
    parser = argparse.ArgumentParser(
        description="The script takes in a file with bin/MAG accessions"
        "(ERZ, GCA, wgsSet) and outputs the accession of the "
        "primary metagenome each came from where possible"
    )
    parser.add_argument(
        "-i",
        "--infile",
        required=True,
        help="Path to the file containing a list of bin/MAG accessions, one per line.",
    )
    parser.add_argument(
        "-o", "--outfile", required=True, help="Name of the outfile to write MAG-assembly pairs."
    )
    parser.add_argument("--debug", action="store_true", help="Print out more information")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_logging(args.debug)
    main(args.infile, args.outfile)
