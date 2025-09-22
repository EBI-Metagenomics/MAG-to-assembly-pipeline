#!/usr/bin/env python
# coding=utf-8

import argparse
import csv
import gzip
import hashlib
import logging
import shutil
from pathlib import Path
from typing import Optional

import requests
from Bio import SeqIO
from botocore.exceptions import ClientError, ParamValidationError
from download_fasta_utils import (
    download_from_ENA_API,
    download_from_ENA_FIRE,
    download_from_ENA_FTP,
)


def main(input_file, output_verified_file, output_invalid_file, download_folder, cleanup):
    validated_pairs = []
    invalid_pairs = []

    with open(input_file, "r") as input_handle:
        reader = csv.reader(input_handle, delimiter="\t")
        for row in reader:
            genome = row[0]
            assemblies = row[1].split(",")
            logging.info(f"Start processing of genome {genome}")
            logging.debug("Verify retrieved assemblies using comparason of contigs' hashes")
            mag_hashes = handle_fasta_processing(genome, download_folder, write_cache=False)
            if not mag_hashes:
                logging.info(f"Failed to download genome {genome} fasta file. Skipping")
                continue

            logging.debug("Genome hashes were computed")
            logging.debug("Start comparing genome hashes to the hashes of the primary assembly")
            confirmed_assemblies = []
            not_confirmed_assemblies = []
            for assembly in assemblies:
                assembly_hashes = handle_fasta_processing(
                    assembly, download_folder, write_cache=True
                )
                if not assembly_hashes:
                    # TODO: think if this case should be considered as an error
                    logging.info(
                        f"For the genome {genome} failed to download primary assembly {assembly} fasta file"
                    )
                    continue
                logging.debug("Assembly hashes were computed")
                if mag_hashes.issubset(assembly_hashes):
                    logging.debug(
                        f"Assembly {assembly} is confirmed to be primary assembly for the genome {genome}"
                    )
                    confirmed_assemblies.append(assembly)
                else:
                    logging.debug(
                        f"Assembly {assembly} is not a primary assembly for the genome {genome}"
                    )
                    not_confirmed_assemblies.append(assembly)
            logging.debug("Comparason finished")

            if confirmed_assemblies:
                logging.info(
                    f"Genome {genome} has been validated to originate from assemblies: {','.join(confirmed_assemblies)}"
                )
                validated_pairs.append([genome, ",".join(confirmed_assemblies)])
            else:
                logging.info(
                    f"Genome {genome} does not have any assemblies with matching contig hashes"
                )
                invalid_pairs.append(
                    [
                        genome,
                        f"the contigs of this genome are not identical to the contigs in the primary assemblies {','.join(assemblies)}",
                    ]
                )

    logging.debug(
        f"Writing validated genome - assembly pairs to the output file {output_verified_file}"
    )
    with open(output_verified_file, "w") as out:
        writer = csv.writer(out, delimiter="\t")
        for line in validated_pairs:
            writer.writerow(line)

    logging.debug(
        f"Writing invalid genome - assembly pairs to the output file {output_invalid_file}"
    )
    with open(output_invalid_file, "w") as out:
        writer = csv.writer(out, delimiter="\t")
        for line in invalid_pairs:
            writer.writerow(line)

    if cleanup and download_folder.exists():
        shutil.rmtree(download_folder)
        logging.debug("Folder with downloaded fasta files is deleted")


def handle_fasta_processing(
    accession: str, download_folder: Path, write_cache=False
) -> Optional[set]:
    """
    Downloads fasta file for the given accession, using different approaches depending on the accession type.
    ERZ analysis accessions are downloaded using ENA FIRE API (fastest and most reliable method, only works for ERZ).
    WGS set accessions are downloaded using ENA FTP links (less reliable method).
    GCA accessions are downloaded using ENA fasta download API (least reliable, but the only good way to download GCA).
    Then computes and returns set of md5 hashes of contigs.
    If the file is already downloaded and cached, it reads hashes from the cache.

    :param accession: ENA accession (e.g. GCA_000001405.15, ERZ1234567)
    :param download_folder: Folder to store downloaded files
    :param write_cache: If True, writes computed hashes to a cache file
    :return: set of md5 hashes of contigs
    """
    if not download_folder.exists():
        download_folder.mkdir(parents=True)
        logging.debug(f"Directory {download_folder} is created")

    try:
        outpath = download_folder / f"{accession}.fa.gz"
        cache_path = download_folder / f"{accession}.hash"
        if (outpath.exists() and outpath.stat().st_size != 0) or (
            cache_path.exists() and cache_path.stat().st_size != 0
        ):
            return compute_hashes(outpath, write_cache=write_cache)

        if accession.startswith("ERZ"):
            # In ENA generated_ftp or submitted_ftp (or both) fields may contain invalid links
            # so we try to download from generated_ftp first, then from submitted_ftp
            try:
                fasta_file = download_from_ENA_FIRE(accession, "generated_ftp", outpath)
                return compute_hashes(fasta_file, write_cache=write_cache)
            except (gzip.BadGzipFile, ClientError, ParamValidationError, ValueError) as e:
                logging.error(
                    f"{accession} Download from link in 'generated_ftp' failed due to: {e}"
                )
                logging.debug('Retry with "submitted_ftp"')
                fasta_file = download_from_ENA_FIRE(accession, "submitted_ftp", outpath)
                return compute_hashes(fasta_file, write_cache=write_cache)

        elif accession.startswith("GCA"):
            fasta_file = download_from_ENA_API(accession, outpath)
            return compute_hashes(fasta_file, write_cache=write_cache)

        else:
            fasta_file = download_from_ENA_FTP(accession, outpath)
            if fasta_file is None:
                raise ValueError("Empty URL or empty file")
            return compute_hashes(fasta_file, write_cache=write_cache)

    except requests.HTTPError as e:
        status_code = e.response.status_code if e.response else "?"
        reason = e.response.reason if e.response else "?"
        logging.error(f"{accession} HTTP Error while downloading: {status_code} - {reason}")
        return None
    # TODO: remove or handle more specific exceptions if needed
    except Exception as e:
        logging.error(f"{accession} Failed to process fasta file due to: {e}")
        return None


def compute_hashes(file_path, write_cache=True, delete_fasta=True, cache_dir=None) -> set:
    """
    Computes md5 hashes of sequences in a fasta file. Creates a cache file with hashes if write_cache is True.
    If the cache file already exists, reads hashes from it.
    :param file_path: Path to the fasta file (can be compressed with .gz)
    :param write_cache: If True, writes computed hashes to a cache file
    :param delete_fasta: If True, deletes the fasta file after processing
    :param cache_dir: Directory to store the cache file, if None uses the same directory as the fasta file
    :return: set of md5 hashes of sequences
    """
    file_path = Path(file_path)
    hashes = set()
    cache_path = file_path.with_suffix("").with_suffix(".hash")
    if cache_dir:
        cache_path = Path(cache_dir) / cache_path.name

    if cache_path.exists():
        with cache_path.open("r") as handle:
            for line in handle:
                hashes.add(line.strip())
        return hashes

    if file_path.suffix == ".gz":
        with gzip.open(file_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                hash_object = hashlib.md5(str(record.seq.upper()).encode())
                hashes.add(hash_object.hexdigest())
    else:
        for record in SeqIO.parse(file_path, "fasta"):
            hash_object = hashlib.md5(str(record.seq.upper()).encode())
            hashes.add(hash_object.hexdigest())

    if not hashes:
        raise ValueError("Fasta file does not contain any records")

    if write_cache:
        with cache_path.open("w") as handle:
            for hash in hashes:
                handle.write(hash + "\n")

    if delete_fasta:
        file_path.unlink()

    return hashes


def parse_args():
    parser = argparse.ArgumentParser(
        description="The script takes in a TSV file with genome - assemblies relationships "
        "(assemblies are separated by commas in the same column), "
        "downloads FASTA files of the assemblies and genomes from ENA, "
        "computes checksums of all contigs, and verifies "
        "if the checksums of the genome contigs match the checksums of the assembly contigs. "
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to the file containing a mapping of genomes to assemblies",
    )
    parser.add_argument(
        "--output_verified",
        required=True,
        type=Path,
        help="Name of the output file to write genome-assembly pairs where all genome contigs match assembly contigs",
    )
    parser.add_argument(
        "--output_invalid",
        required=True,
        type=Path,
        help="Name of the output file to write genome-assembly pairs that failed contig validation",
    )
    parser.add_argument(
        "--errors",
        required=True,
        type=Path,
        help="Name of the output file to log ENA-ralated errors during fasta download",
    )
    parser.add_argument(
        "--download-folder",
        required=False,
        type=Path,
        help="Folder to store downloaded files. By default: fasta_downloads",
        default="downloaded_fastas",
    )
    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Remove downloaded cache of checksums and download folder after execution",
    )
    parser.add_argument("--debug", action="store_true", help="Increase logging verbosity")
    return parser.parse_args()


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
    for noisy_lib in ["requests", "boto3", "botocore", "urllib", "urllib3", "s3transfer"]:
        logging.getLogger(noisy_lib).setLevel(logging.WARNING)


if __name__ == "__main__":
    args = parse_args()
    setup_logging(args.debug, args.errors)
    main(args.input, args.output_verified, args.output_invalid, args.download_folder, args.cleanup)
