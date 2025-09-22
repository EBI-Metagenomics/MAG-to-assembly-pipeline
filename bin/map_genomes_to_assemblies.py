#!/usr/bin/env python
# coding=utf-8

import argparse
import csv
import logging
import re
import time
from functools import wraps
from typing import Optional

import requests
import xmltodict

# TODO look for primary assemblies even if bin sample is bio sample?


def main(input_file, output_file, no_assembly_file):
    """
    Main function to map bin/MAG genome accessions to primary metagenome assemblies.
    It goes through the following steps for each genome accession:
    1. Determine the sample accession associated with the genome accession.
    2. Retrieve the "derived from" sample accessions from the genome's sample metadata
    3. Find all primary metagenome assemblies linked to the "derived from" samples.
    4. If multiple assemblies are found, filter them based on run accessions if available.

    :param input_file: Path to the input file containing genome accessions, one per line.
    :param output_file: Path to the output file to write genome-assemblies mapping.
    :param no_assembly_file: Path to the output file to log genomes with no assemblies found.
    """
    genome_assembly_pairs = []
    no_assembly_genomes = []

    with open(input_file, "r") as file_in:
        reader = csv.reader(file_in)
        for row in reader:
            genome_accession = row[0].strip()
            if genome_accession[:3] not in ["ERZ", "GCA"]:
                genome_accession = genome_accession.rstrip("0")  # CAMPAA010000000 -> CAMPAA01

            logging.debug(f"Start processing of genome {genome_accession}")
            logging.debug(f"Find sample accession corresponding to the genome {genome_accession}")
            try:
                genome_sample = find_genome_sample_in_ena(genome_accession)
            except (requests.exceptions.RequestException, KeyError, ValueError) as e:
                log_skipped_genome(
                    genome_accession,
                    f"Failed to fetch sample accession for the genome {genome_accession} due to {e}",
                    no_assembly_genomes,
                )
                continue
            logging.debug(
                f"Successful. Genome {genome_accession} sample accession is {genome_sample}"
            )

            logging.debug(
                f"Find 'derived from' sample corresponding to the genome {genome_accession}"
            )
            try:
                derived_from_samples, derived_from_runs = find_derived_from_sample_in_ena(
                    genome_sample
                )
            except requests.exceptions.RequestException as e:
                log_skipped_genome(
                    genome_accession,
                    f"Unable to get genome sample XML for {genome_sample} due to: {e}",
                    no_assembly_genomes,
                )
                continue
            except (KeyError, ValueError) as e:
                log_skipped_genome(
                    genome_accession,
                    f"Failed to parse XML for {genome_sample} due to: {e}",
                    no_assembly_genomes,
                )
                continue
            if not derived_from_samples:
                logging.debug("Failed to find sample accessions in the 'derived from' attribute")
                logging.debug(
                    "Attempting to collect sample accessions from the runs' metadata using ENA API"
                )
                try:
                    derived_from_samples = collect_samples_from_runs_metadata(derived_from_runs)
                    if not derived_from_samples:
                        raise ValueError("No sample accessions found in run metadata")
                except (requests.exceptions.RequestException, KeyError, ValueError) as e:
                    log_skipped_genome(
                        genome_accession,
                        f"Unable to find sample accessions from runs for {genome_sample} due to {e}",
                        no_assembly_genomes,
                    )
                    continue
            logging.debug(
                f"The following sample accessions were found: {comma_separate(derived_from_samples)}"
            )
            # TODO sometimes for ERZ genomes runs references are included in the analysis XML
            # look for runs in the assembly XML if genome is an assembly
            if not derived_from_runs:
                logging.debug(
                    f"Runs for the genome {genome_accession} were not found. "
                    f"Comparson of run accessions for the genome and primary assemblies will be skipped"
                )

            logging.debug(
                f"Find all primary metagenomic assemblies linked to the samples {comma_separate(derived_from_samples)}"
            )
            primary_assemblies, assembly2runs = get_primary_assemblies_from_sample(
                derived_from_samples
            )
            if not primary_assemblies:  # cases when primary assembly was not uploaded to ENA
                log_skipped_genome(
                    genome_accession,
                    f"There are no assemblies for 'derived from' samples {comma_separate(derived_from_samples)}",
                    no_assembly_genomes,
                )
                continue
            logging.debug(
                f"Successful. The following primary assemblies were found {comma_separate(primary_assemblies)}"
            )

            if len(primary_assemblies) > 1 and derived_from_runs:
                logging.debug(
                    "Attempt to decrease the list of assemblies by filtering assemblies "
                    "generated from the runs other than genome's runs"
                )
                for assembly, runs in assembly2runs.items():
                    if runs and runs != set(derived_from_runs):
                        primary_assemblies.remove(assembly)
                if not primary_assemblies:
                    log_skipped_genome(
                        genome_accession,
                        f"There are no assemblies with similar runs for 'derived from' sample {comma_separate(derived_from_samples)}",
                        no_assembly_genomes,
                    )
                    continue
                logging.debug(f"Updated list of assemblies: {comma_separate(primary_assemblies)}")

            if primary_assemblies:
                genome_assembly_pairs.append([genome_accession, comma_separate(primary_assemblies)])
            else:
                log_skipped_genome(
                    genome_accession,
                    f"No primary assemblies found for 'derived from' samples {comma_separate(derived_from_samples)}",
                    no_assembly_genomes,
                )

    logging.debug(f"Write list of identified assemblies to the file {output_file}")
    with open(output_file, "w") as file_out:
        writer = csv.writer(file_out, delimiter="\t")
        for line in genome_assembly_pairs:
            writer.writerow(line)

    logging.debug(f"Write list of genomes with no assembly found to the file {no_assembly_file}")
    with open(no_assembly_file, "w") as file_out:
        writer = csv.writer(file_out, delimiter="\t")
        for line in no_assembly_genomes:
            writer.writerow(line)


def log_skipped_genome(genome_accession: str, message: str, log_list: list) -> None:
    logging.debug(message)
    log_list.append([genome_accession, message])
    logging.debug("Skipping")


def find_genome_sample_in_ena(genome_accession: str) -> str:
    """
    Given a genome accession (ERZ, GCA, wgsSet), retrieve its metadata from ENA portal
    and parse it to get the associated sample accession.
    """
    if genome_accession[:3] in ["ERZ", "GCA"]:
        logging.debug(
            f"{genome_accession} is an ENA analysis or genome accession, "
            f"retrieving metadata in XML from ENA portal"
        )
        genome_ena_data = load_data(genome_accession, endpoint_type="xml")
        accession_type = "ANALYSIS" if genome_accession.startswith("ERZ") else "ASSEMBLY"
        sample_ref = genome_ena_data[f"{accession_type}_SET"][f"{accession_type}"]["SAMPLE_REF"]
        return sample_ref["IDENTIFIERS"]["PRIMARY_ID"]
    else:
        logging.debug(
            f"{genome_accession} is an ENA WGS set accession, retrieving summary from ENA portal"
        )
        genome_ena_data = load_data(genome_accession, endpoint_type="summary")
        return genome_ena_data["summaries"][0]["sample"]


def find_derived_from_sample_in_ena(genome_sample: str) -> tuple:
    """
    Given a sample accession, retrieve its metadata in XML from ENA portal and parse the "derived from" attribute.
    Return a tuple of two elements:
    1. A list of sample accessions from which the given sample was derived (if available).
    2. A list of run accessions from which the given sample was derived (if available).
    If neither is found, raise a ValueError.
    """
    logging.debug(f"Retrieving metadata in XML for accession {genome_sample} from ENA portal")
    sample_data = load_data(genome_sample, endpoint_type="xml")
    sample_attributes = sample_data["SAMPLE_SET"]["SAMPLE"]["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]
    logging.debug("Parsing sample attributes in XML metadata")
    derived_from_samples = []
    derived_from_runs = []

    for attribute in sample_attributes:
        # TODO can this attribute have a different name?
        if all([x in attribute["TAG"] for x in ["derived", "from"]]):
            derived_from_samples.extend(
                re.findall(r"SAM[A-Z]+\d+|ERS\d+|SRS\d+|DRS\d+", attribute["VALUE"])
            )
            derived_from_runs.extend(re.findall(r"ERR\d+|SRR\d+|DRR\d+", attribute["VALUE"]))
            break
    if not derived_from_runs and not derived_from_samples:
        raise ValueError(f"No 'derived from' attribute in sample XML {genome_sample}")

    return derived_from_samples, derived_from_runs


def get_primary_assemblies_from_sample(sample_accessions: list) -> tuple:
    """
    Given a list of sample accessions, find all primary metagenome assemblies linked to these
    samples using ENA search API.
    Return a tuple of two elements:
    1. A list of primary assembly accessions.
    2. A dictionary mapping each assembly accession to its corresponding run accession(s).
    """
    primary_assemblies = []
    assembly2runs = {}
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
            run_accessions = set(row["run_accession"].split(","))
            primary_assemblies.append(assembly_accession)
            assembly2runs[assembly_accession] = run_accessions

    return sorted(primary_assemblies), assembly2runs


def genbank_to_ena_wgsset_accession(genome_accession: str) -> Optional[str]:
    """
    Convert a GenBank WGS accession (GCA) to the corresponding ENA WGS set accession (wgsSet).
    """
    try:
        summary_data = load_data(genome_accession, endpoint_type="summary")
        ena_accession = summary_data["summaries"][0]["wgsSet"]
        return ena_accession
    except (requests.exceptions.RequestException, KeyError, ValueError) as e:
        logging.info(
            f"Unable to convert NCBI accession {genome_accession} to wgsSet accession due to {e}"
        )
        return None


def load_data(accession: str, endpoint_type: str) -> dict:
    """
    Retrieve metadata for a given ENA accession from XML or summary endpoints.

    :param accession: ENA accession
    :param endpoint_type: 'xml' or 'summary' indicating the API endpoint to use
    :return: Metadata as a dictionary
    """
    url = f"https://www.ebi.ac.uk/ena/browser/api/{endpoint_type}/{accession}"
    try:
        response = run_browser_request(url)
        if endpoint_type == "xml":
            return xmltodict.parse(response.content)
        elif endpoint_type == "summary":
            return response.json()
        else:
            raise ValueError(
                f"Unsupported endpoint_type '{endpoint_type}'. Must be 'xml' or 'summary'."
            )
    except requests.exceptions.RequestException as e:
        logging.error(f"{accession} unable to request content from {url}: {e}")
        raise


def collect_samples_from_runs_metadata(runs: list) -> list:
    """
    Given a list of run accessions, retrieve their metadata in XML from ENA portal
    and parse it to get sample accessions.
    """
    samples = []
    for run in runs:
        run_data = load_data(run, endpoint_type="xml")
        for link in run_data["RUN_SET"]["RUN"]["RUN_LINKS"]["RUN_LINK"]:
            if link["XREF_LINK"]["DB"] == "ENA-SAMPLE":
                sample = link["XREF_LINK"]["ID"]
                samples.append(sample)
                break
    return samples


def custom_retry(max_retries=3, delay=10, backoff=1.5):
    """Retry decorator with special handling for 404 errors.
    - 404: retry only once
    - Other errors: retry up to max_retries with backoff
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            attempt = 0
            current_delay = delay
            while True:
                try:
                    return func(*args, **kwargs)
                except requests.exceptions.HTTPError as e:
                    status = e.response.status_code
                    attempt += 1
                    if status == 404:
                        if attempt > 1:
                            logging.error(f"404 Not Found, giving up after {attempt} attempt(s).")
                            raise
                        logging.warning(
                            f"404 Not Found, retrying once in {current_delay} seconds..."
                        )
                        time.sleep(current_delay)
                        continue
                    else:
                        if attempt > max_retries:
                            logging.error(
                                f"HTTP error {status}, giving up after {attempt} attempts."
                            )
                            raise
                        logging.warning(
                            f"HTTP error {status}, retrying in {current_delay:.1f} seconds..."
                        )
                        time.sleep(current_delay)
                        current_delay *= backoff
                except requests.exceptions.RequestException as e:
                    attempt += 1
                    if attempt > max_retries:
                        logging.error(f"Request failed, giving up after {attempt} attempts: {e}")
                        raise
                    logging.warning(
                        f"Request failed: {e}, retrying in {current_delay:.1f} seconds..."
                    )
                    time.sleep(current_delay)
                    current_delay *= backoff

        return wrapper

    return decorator


@custom_retry(max_retries=5, delay=15, backoff=1.5)
def run_browser_request(url):
    """Run a GET request to the ENA browser API."""
    response = requests.get(url)
    response.raise_for_status()
    return response


@custom_retry(max_retries=5, delay=15, backoff=1.5)
def run_request(query, api_endpoint):
    """Run a GET request to the ENA search API with given query parameters."""
    response = requests.get(api_endpoint, params=query)
    response.raise_for_status()
    return response


def setup_logging(debug=False, error_logfile="ena_related_errors.log"):
    """
    Set up logging configuration.
    :param debug: If True, set logging level to DEBUG, otherwise to INFO
    :param error_logfile: Path to the file where errors appearing during fasta download will be logged
    """
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


def comma_separate(items: list) -> str:
    """Convert a list of items into a comma-separated string."""
    return ",".join(items)


def parse_args():
    parser = argparse.ArgumentParser(
        description="The script takes in a file with bin/MAG accessions"
        "(ERZ, GCA, wgsSet) and outputs the accession of the "
        "primary metagenome each came from where possible"
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the file containing a list of bin/MAG genome accessions, one per line.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Name of the output_file to write genome-assembly pairs.",
    )
    parser.add_argument(
        "--no_assembly_found",
        required=True,
        help="Name of the output file to write genome accessions for which no primary assembly was found.",
    )
    parser.add_argument(
        "-e",
        "--errors",
        required=True,
        help="Name of the output file to log ENA-related errors.",
    )
    parser.add_argument("--debug", action="store_true", help="Print out more information")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_logging(args.debug, args.errors)
    main(args.input, args.output, args.no_assembly_found)
