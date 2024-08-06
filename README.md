# MAG-to-Assembly Matching Pipeline

This repository contains an in-progress Nextflow pipeline designed to match Metagenome-Assembled Genome (MAG) accessions to their corresponding primary metagenome assemblies. The pipeline uses metadata links retrieved through the ENA API and verifies matches by comparing contig checksums between MAGs and assemblies.

The pipeline consists of three consecutive processes:

1. **Creation of Input for the Main Script**

   Download lists of accessions from MGnify catalogs and ENA, merge them, remove redundancy, and exclude accessions that have been processed in previous runs of the pipeline.

2. **MAG-to-Assembly Linking**

   This process utilizes the ENA API to retrieve metadata related to the provided MAG accessions and download FASTA files for both MAGs and assemblies. These files are necessary to compute checksum values of the contigs.

3. **Postprocessing**

   Merge results into a single table, add `Species_rep` and `Action` columns, and create updated `processed_accessions_*.tsv` and `mag_to_assembly_links_*.tsv` files.

## Requirements

- Singularity
- Nextflow

## Input Arguments
 ### Required 
- `--previous_table`: Earlier `mag_to_assembly_links_*.tsv` file generated by the pipeline in a previous run. An example is `data/test/mag_to_assembly_links.tsv`
- `--processed_acc`: List of the accessions `processed_accessions_*.tsv` that have been processed in previous runs of the pipeline. An example is `data/test/already_processed_accessions.tsv`

__If these arguments are not specified, the pipeline will run inference of MAG-to-assembly links for all MAGs and bins deposited in the ENA and MGnify catalogues.__

### Optional 
- `output_path`: Where to save output files. Default is `$PWD`

## Usage

You can run the pipeline on SLURM with the following command (use absolute paths):

```bash
nextflow run -profile codon_slurm main.nf --previous_table /path/to/file.tsv --processed_acc /path/to/file.tsv
```

To launch a run on a test dataset `data/test/`, use:

```bash
nextflow run -profile codon_slurm main.nf --test=true
```

