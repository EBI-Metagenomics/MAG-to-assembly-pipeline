#!/usr/bin/env nextflow

/// python3 ~/various-data-processing-scripts/link_mags_to_assemblies/link_MAG_to_primary_metagenome_assembly.py -i gut_genomes.mags.tsv -o gut.confirmed.tsv -f gut.fails.tsv -p gut.putative.tsv --download-folder fasta_files

params.download_folder = 'fasta_files'
params.prefix = 'sample'
params.input_file = ""

params.output_file = "${params.prefix}.confirmed.tsv"
params.putative_file = "${params.prefix}.putative.tsv"
params.failed_file = "${params.prefix}.failed.tsv"

process RETROFIT {

    input:
    path input_file

    output:
    path "*.confirmed.tsv"

    script:
    """
    link_MAG_to_primary_metagenome_assembly.py -i ${input_file} -o ${params.prefix}.confirmed.tsv -p ${params.prefix}.putative.tsv -f ${params.prefix}.failed.tsv --download-folder ${params.download_folder}
    """
}

workflow {
    accessions_list_ch = Channel.fromPath(params.input_file, checkIfExists: true)

    // mag_list_ch = DOWNLOAD(accessions_list_ch)
    mag_assembly_pairs_ch = RETROFIT(accessions_list_ch)

    // mag_assembly_pairs_ch.output.putative_file
    // POSTPROCESSING(mag_assembly_pairs_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}
