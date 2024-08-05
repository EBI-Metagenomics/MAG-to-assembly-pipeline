#!/usr/bin/env nextflow

include { DOWNLOAD_INPUT                } from '../modules/download_input_acc'
include { FIND_PRIMARY_ASSEMBLY         } from '../modules/find_assembly'
include { POSTPROCESSING                } from '../modules/finalise_output'

workflow MAG_ASSEMBLY_LINKING_PIPELINE {

    // Check if test mode is enabled and use test values if true
    if (params.test) {
        DOWNLOAD_INPUT(params.test_processed_acc, params.input_accessions, params.gut_mapping, params.catalogue_metadata)
        accessions_list_ch = Channel.fromPath(params.test_input)
    } else {
        DOWNLOAD_INPUT(params.processed_acc, params.input_accessions, params.gut_mapping, params.catalogue_metadata)
        accessions_list_ch = DOWNLOAD_INPUT.output.input_accessions
    } 

    // calculatePortionSize(accessions_list_ch, params.n_tasks)
    accessions_portions_ch = accessions_list_ch.splitText(by: params.portion_size, file: "portion")
    
    FIND_PRIMARY_ASSEMBLY(accessions_portions_ch)
    mag_assembly_pairs_ch = FIND_PRIMARY_ASSEMBLY.output.mag_assembly_pairs
    not_linked_mags_ch = FIND_PRIMARY_ASSEMBLY.output.not_linked_mags


    if (params.test) {
        previous_table_ch = Channel.fromPath(params.test_previous_table)
        processed_acc_ch = Channel.fromPath(params.test_processed_acc)
    } else {
        previous_table_ch = Channel.fromPath(params.previous_table)
        processed_acc_ch = Channel.fromPath(params.processed_acc)
    }

    final_ch = POSTPROCESSING(mag_assembly_pairs_ch.collect(), not_linked_mags_ch.collect(), DOWNLOAD_INPUT.output.metadata, processed_acc_ch, previous_table_ch)
}

