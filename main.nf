#!/usr/bin/env nextflow

// DOWNLOAD input args
params.processed_acc = "$projectDir/data/previously_processed_acc.test.tsv"
params.gut_mapping = "$projectDir/data/known_and_found_gut_genomes.tsv"
params.input_accessions = 'input_accessions.tsv'
params.catalogue_metadata = 'all_catalog_metadata.tsv'
// params.input_accessions = ""


// RETROFIT input args
params.portion_size = 10 

// POSTPROCESSING input args
params.previous_table = "$projectDir/data/last_version_table.tsv"
// params.


process DOWNLOAD {
    input:
    path processed_acc
    val input_accessions
    path gut_mapping
    val catalogue_metadata
    

    output:
    path "${input_accessions}", emit: input_accessions
    path catalogue_metadata, emit: metadata

    script:
    """
    preprocessing.py --processed-acc ${processed_acc} --output-accessions ${input_accessions} --gut-mapping ${gut_mapping} --catalogue-metadata ${catalogue_metadata}
    """
}


process RETROFIT {
    input:
    path input

    output:
    path "*.confirmed.tsv"

    script:
    """
    link_MAG_to_primary_metagenome_assembly.py -i ${input} -o ${input}.confirmed.tsv -p ${input}.putative.tsv -f ${input}.failed.tsv --download-folder fastas
    """
}

process POSTPROCESSING {
    input:
    path input
    path catalogue_metadata

    output:
    path 'RETROFIT.tsv'

    script:
    """
    postprocessing.py --previous-table ${params.previous_table} --catalogue-metadata ${catalogue_metadata} ${input}
    """
}
workflow {

    DOWNLOAD(params.processed_acc, params.input_accessions, params.gut_mapping, params.catalogue_metadata)
    // input_accessions_ch.view()

    accessions_list_ch = DOWNLOAD.output.input_accessions.splitText(by: params.portion_size, file: "portion")
    // accessions_list_ch.view()

    mag_assembly_pairs_ch = RETROFIT(accessions_list_ch)
    mag_assembly_pairs_ch.view()

    final_ch = POSTPROCESSING(mag_assembly_pairs_ch.collectFile(), DOWNLOAD.output.metadata)
    final_ch.view()
}

// workflow.onComplete {
//     log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
// }
