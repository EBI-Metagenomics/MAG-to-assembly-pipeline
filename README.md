# MAG-to-Assembly Linking Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

This Nextflow pipeline is designed to map Metagenome-Assembled Genome (MAG) accessions to their corresponding primary metagenome assemblies. The pipeline retrieves metadata links through the [ENA Portal API and Browser API](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access.html), than verifies matches by comparing contig checksums between MAGs and assemblies.

## Pipeline Overview

The pipeline performs the following main steps:

1. **Creation of the list of input genome accessions**

   Gather accessions from MGnify catalogs and ENA, merge them, remove redundancy, and exclude accessions that have been processed in previous runs of the pipeline. This step is skipped if running on user-provided input accessions.

2. **Mapping of genomes to assemblies**

   This process utilizes the ENA API to retrieve metadata of the input genome accessions and, using connections between genomes, samples, and assemblies in the ENA data model, identifies accessions of putative primary assemblies.

3. **Validation of found assemblies using contig checksums**

   Genome-assembly pair is only considered valid if all genome's contigs are a subset of assembly's contigs. This is validated through comparison of hash sets of contig sequences for each assembly and genome. To compute contig checksum values, FASTA files for all genomes and assemblies are downloaded.

4. **Formatting of output files**

   MAG-assembly pairs are merged into a single table, `Species_rep` column is added, and updated [`processed_accessions_*.tsv`](workflows/tests/data/processed_accessions.tsv) and [`mag_to_assembly_mapping_*.tsv`](workflows/tests/data/mag_to_assembly_links.tsv) files are created.

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=23.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)). You can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles).

3. Clone the repository:

   ```bash
   git clone https://github.com/EBI-Metagenomics/MAG-to-assembly-pipeline.git
   cd MAG-to-assembly-pipeline
   ```

4. Run the pipeline with test data on a MacOS machine using Docker containers:

   ```bash
   nextflow run main.nf \
     --input_accessions workflows/tests/data/input_accessions.tsv \
     --catalogues_metadata workflows/tests/data/all_catalog_metadata.tsv \
     --merge_with_results workflows/tests/data/mag_to_assembly_links.tsv \
     -profile docker,test,arm
   ```

## Usage

By default, the pipeline does not require any input from the user.

### Optional Input Parameters

- `--accessions_list`: Path to TSV file containing input genome accessions (one per line). Use this to process a custom list of accessions instead of those collected automatically in the first step of the pipeline.
- `--catalogues_metadata`: Path to TSV file containing genomes' metadata. Only used if provided with `--accessions_list`.
- `--merge_with_results`: Path to existing genome-assembly mapping file to merge with new results (see [example](workflows/tests/data/mag_to_assembly_links.tsv)).
- `--skip_accessions`: Path to TSV file containing accessions (one per line) that have been processed in previous runs and should be excluded from processing. Enabled when `--accessions_list` is used.

### Optional Output Parameters

- `--output_path`: Output directory where results will be saved (default: `./results`)

### Optional Execution Parameters

- `--batch_size`: Size of batches into which the list of input accessions is divided for parallel processing (default: `250`)
- `--debug`: Enable debug mode to generate additional logging information (default: `false`)
- `--cleanup`: Remove contig hashes cache after processing to save disk space (default: `true`)

### Example Commands

#### Run with default settings (process all ENA/MGnify accessions)

```bash
nextflow run main.nf -profile docker
```

#### Run with custom accession list

```bash
nextflow run main.nf \
  --accessions_list input_accessions.tsv \
  --catalogues_metadata metadata.tsv \
  --output_path results \
  -profile docker
```

#### Run with previous results to merge

```bash
nextflow run main.nf \
  --merge_with_results previous_mapping.tsv \
  --skip_accessions processed_accessions.tsv \
  -profile docker
```

## Output

The pipeline generates the following outputs in the specified output directory:

- `processed_accessions_YYYY-MM-DD_HHhMMm.tsv`: Timestamped file containing processed accession results
- `mag_to_assembly_mapping_YYYY-MM-DD_HHhMMm.tsv`: Timestamped file containing mapping of genomes to assemblies
- `ena_related_errors/*.err`: Error logs from processing steps listing genomes that failed processing due to ENA-related issues.

## Configuration Profiles

The pipeline comes with several configuration profiles:

- `local`: For running on local machines with limited resources
- `test`: For running in the test mode
- `arm`: For running on MacOS machines
- `codon_slurm`: For running on SLURM clusters (specifically configured for Codon)
- `docker`: For running with Docker containers
- `singularity`: For running with Singularity containers
