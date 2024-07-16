#!/usr/bin/env nextflow

/// python3 ~/various-data-processing-scripts/link_mags_to_assemblies/link_MAG_to_primary_metagenome_assembly.py -i gut_genomes.mags.tsv -o gut.confirmed.tsv -f gut.fails.tsv -p gut.putative.tsv --download-folder fasta_files

params.working_folder = ""
params.download_folder = ""
params.prefix = ""
params.output_file = "${params.prefix}.confirmed.tsv"
params.putative_file = "${params.prefix}.putative.tsv"
params.failed_file = "${params.prefix}.failed.tsv"

process DOWNLOAD {
    input:
    path previous_retrofit

    output:
    path 'mags_list.txt'

    script:
    """
    python download_script.py --input ${previous_retrofit} --output mags_list.txt
    """
}

process RETROFIT {
    input:
    path 'mags_list.txt'

    output:
    path 'mag_assembly_pairs.txt'

    script:
    """
    python retrofit_script.py --input mags_list.txt --output mag_assembly_pairs.txt
    """
}

process POSTPROCESSING {
    input:
    path 'mag_assembly_pairs.txt'

    output:
    path 'final_retrofit.txt'

    script:
    """
    python postprocessing_script.py --input mag_assembly_pairs.txt --output final_retrofit.txt
    """
}


workflow {
    previous_retrofit_ch = Channel.fromPath(params.previous_retrofit, checkIfExists: true)

    mag_list_ch = DOWNLOAD(previous_retrofit_ch)
    mag_assembly_pairs_ch = RETROFIT(mag_list_ch)
    POSTPROCESSING(mag_assembly_pairs_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}
