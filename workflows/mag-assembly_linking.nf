#!/usr/bin/env nextflow

include { DOWNLOAD_INPUT                } from '../modules/download_input'
include { FIND_PRIMARY_ASSEMBLY         } from '../modules/find_assembly'
include { FINALISE_OUTPUT               } from '../modules/finalise_output'

workflow MAG_ASSEMBLY_LINKING_PIPELINE {

    processed_acc_ch = params.processed_acc ? Channel.fromPath(params.processed_acc) : []

    // If custom input accessions are provided, use them instead of the downloaded accessions
    if (params.external_input) {
        accessions_list_ch = Channel.fromPath(params.external_input)
        metadata_ch = params.external_metadata ? Channel.fromPath(params.external_metadata) : []

    } else {
        DOWNLOAD_INPUT(processed_acc_ch, params.input_accessions, params.gut_mapping, params.catalogue_metadata)
        metadata_ch = DOWNLOAD_INPUT.output.metadata
        accessions_list_ch = DOWNLOAD_INPUT.output.input_accessions
    }

    // Input accessions are splitted to process them faster in parallel tasks
    accessions_portions_ch = accessions_list_ch.splitText(by: params.portion_size, file: "portion")
    FIND_PRIMARY_ASSEMBLY(accessions_portions_ch)

    mag_assembly_pairs_ch = FIND_PRIMARY_ASSEMBLY.output.mag_assembly_pairs
    not_linked_mags_ch = FIND_PRIMARY_ASSEMBLY.output.not_linked_mags
    previous_table_ch = params.previous_table ? Channel.fromPath(params.previous_table) : []
    final_ch = FINALISE_OUTPUT(mag_assembly_pairs_ch.collect(), not_linked_mags_ch.collect(), metadata_ch, processed_acc_ch, previous_table_ch)
}

