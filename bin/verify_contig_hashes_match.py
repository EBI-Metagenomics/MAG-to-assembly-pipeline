import csv
import gzip
import hashlib
import logging
import os
import shutil

import requests
from Bio import SeqIO
from botocore.exceptions import ClientError, ParamValidationError
from download_fasta_utils import (
    download_from_ENA_API,
    download_from_ENA_FIRE,
    download_from_ENA_FTP,
)


def main(input_file, download_folder, cleanup, out_confirmed, out_putative, out_fails):

    reader = csv.DictReader(input_file, delimiter="\t")
    for row in reader:
        bin = row["Genome_accession"]
        assemblies = row["Assembly_accession"]
        logging.debug("Verify retrieved assemblies using comparason of contigs' hashes")
        mag_hashes = handle_fasta_processing(bin, download_folder)
        if not mag_hashes:
            logging.info(f"Failed to download MAG {bin} fasta file. Skipping")
            print(bin, f"Failed to download fasta file for bin {bin}.", sep="\t", file=out_fails)
            continue
        logging.debug("Bin hashes were computed")
        logging.debug("Start comparing bin hashes to hashes of the primary assembly")
        confirmed_assemblies = []
        not_confirmed_assemblies = []
        for assembly in assemblies:
            assembly_hashes = handle_fasta_processing(assembly, download_folder)
            if not assembly_hashes:
                logging.info(
                    f"For the MAG {bin} failed to download primary assembly {assembly} fasta file."
                )
                continue
            logging.debug("Assembly hashes were computed")
            if mag_hashes.issubset(
                assembly_hashes
            ):  # TODO modify to avoid matching empty file hashes
                logging.debug(
                    f"Assembly {assembly} is confirmed to be primary assembly for the MAG/bin {bin}"
                )
                confirmed_assemblies.append(assembly)
            else:
                logging.debug(
                    f"Assembly {assembly} is not a primary assembly for the MAG/bin {bin}"
                )
                not_confirmed_assemblies.append(assembly)
        logging.debug("Comparason finished")

        # Write a line to the output TSV file
        # Columns are      Genome_acc    Sample      Derived_from_sample     Derived_from_assembly
        logging.debug("Writing results to the output file")
        if confirmed_assemblies:
            print(bin, ",".join(confirmed_assemblies), sep="\t", file=out_confirmed)
        elif not_confirmed_assemblies:
            print(bin, ",".join(not_confirmed_assemblies), sep="\t", file=out_putative)

        if cleanup and os.path.exists(download_folder):
            shutil.rmtree(download_folder)
            logging.debug("Folder with downloaded files is deleted")


def handle_fasta_processing(accession, download_folder):
    try:
        outpath = os.path.join(download_folder, f"{accession}.fa.gz")
        cache_path = os.path.join(download_folder, f"{accession}.fa.hash")
        if (os.path.exists(outpath) and os.path.getsize(outpath) != 0) or (
            os.path.exists(cache_path) and os.path.getsize(cache_path) != 0
        ):
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
                logging.error(
                    f"{accession} Download from link in 'generated_ftp' failed due to: {e}"
                )
                logging.debug('Retry with "submitted_ftp"')
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

    if not hashes:
        raise ValueError("Fasta file does not contain any records")

    if write_cache:
        with open(cache_path, "w") as handle:
            for hash in hashes:
                handle.write(hash + "\n")

    if delete_fasta:
        os.remove(file_path)

    return hashes
