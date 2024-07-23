#!/usr/bin/env nextflow

include { DOWNLOAD_INPUT                } from '../modules/download_mag_acc'
include { FIND_PRIMARY_ASSEMBLY         } from '../modules/find_assembly'
include { POSTPROCESSING                } from '../modules/finalise_output'

workflow MAG_ASSEMBLY_LINKING_PIPELINE {

    DOWNLOAD_INPUT(params.processed_acc, params.input_accessions, params.gut_mapping, params.catalogue_metadata)
    // input_accessions_ch.view()

    if (params.test) {
        accessions_list_ch = Channel.fromPath(params.test_input)
    } else {
        accessions_list_ch = DOWNLOAD_INPUT.output.input_accessions
    } 

    accessions_portions_ch = accessions_list_ch.splitText(by: params.portion_size, file: "portion")
    // accessions_list_ch.view()
    
    mag_assembly_pairs_ch = FIND_PRIMARY_ASSEMBLY(accessions_portions_ch)
    mag_assembly_pairs_ch.view()

    final_ch = POSTPROCESSING(mag_assembly_pairs_ch.collect(), DOWNLOAD_INPUT.output.metadata)
    final_ch.view()
}

