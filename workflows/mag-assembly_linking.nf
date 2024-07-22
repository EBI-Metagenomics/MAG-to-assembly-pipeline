#!/usr/bin/env nextflow

// DOWNLOAD_INPUT input args
params.processed_acc = "$projectDir/data/previously_processed_acc.test.tsv"
params.gut_mapping = "$projectDir/data/known_and_found_gut_genomes.tsv"
params.input_accessions = 'input_accessions.tsv'
params.catalogue_metadata = 'all_catalog_metadata.tsv'
// params.input_accessions = ""


// FIND_PRIMARY_ASSEMBLY input args
params.portion_size = 10 

// POSTPROCESSING input args
params.previous_table = "$projectDir/data/last_version_table.tsv"
// params.

include { DOWNLOAD_INPUT                } from '../modules/download_mag_acc'
include { FIND_PRIMARY_ASSEMBLY         } from '../modules/find_assembly'
include { POSTPROCESSING                } from '../modules/finalise_output'

workflow MAG_ASSEMBLY_LINKING_PIPELINE {

    DOWNLOAD_INPUT(params.processed_acc, params.input_accessions, params.gut_mapping, params.catalogue_metadata)
    // input_accessions_ch.view()

    accessions_list_ch = DOWNLOAD_INPUT.output.input_accessions.splitText(by: params.portion_size, file: "portion")
    // accessions_list_ch.view()

    mag_assembly_pairs_ch = FIND_PRIMARY_ASSEMBLY(accessions_list_ch)
    mag_assembly_pairs_ch.view()

    final_ch = POSTPROCESSING(mag_assembly_pairs_ch.collect(), DOWNLOAD_INPUT.output.metadata)
    final_ch.view()
}

