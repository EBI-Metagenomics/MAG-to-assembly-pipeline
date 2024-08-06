#!/usr/bin/env nextflow

include { DOWNLOAD_INPUT                } from '../modules/download_input_acc'
include { FIND_PRIMARY_ASSEMBLY         } from '../modules/find_assembly'
include { POSTPROCESSING                } from '../modules/finalise_output'

workflow MAG_ASSEMBLY_LINKING_PIPELINE {

    // Workaround with a empty file channel is required to simulate an optional argument
    processed_acc_ch = params.processed_acc ? Channel.fromPath(params.processed_acc) : Channel.fromPath(params.empty_file)
    DOWNLOAD_INPUT(processed_acc_ch, params.input_accessions, params.gut_mapping, params.catalogue_metadata)

    // If test mode is enabled, instead of the downloaded accessions use the small test file
    if (params.test) {
        accessions_list_ch = Channel.fromPath(params.test_input)
    } else {
        accessions_list_ch = DOWNLOAD_INPUT.output.input_accessions
    }

    // Input accessions are splitted to processes them faster in parallel tasks
    accessions_portions_ch = accessions_list_ch.splitText(by: params.portion_size, file: "portion")
    FIND_PRIMARY_ASSEMBLY(accessions_portions_ch)

    mag_assembly_pairs_ch = FIND_PRIMARY_ASSEMBLY.output.mag_assembly_pairs
    not_linked_mags_ch = FIND_PRIMARY_ASSEMBLY.output.not_linked_mags
    // Workaround with a empty file channel is required to simulate an optional argument
    previous_table_ch = params.previous_table ? Channel.fromPath(params.previous_table) : Channel.fromPath(params.empty_file)
    final_ch = POSTPROCESSING(mag_assembly_pairs_ch.collect(), not_linked_mags_ch.collect(), DOWNLOAD_INPUT.output.metadata, processed_acc_ch, previous_table_ch)
}

