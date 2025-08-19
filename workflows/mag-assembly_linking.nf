#!/usr/bin/env nextflow

include { COLLECT_INPUT_ACCESSIONS   } from '../modules/local/collect_input_accessions/main'
include { MAP_GENOMES_TO_ASSEMBLIES  } from '../modules/local/map_genomes_to_assemblies/main'
include { VERIFY_CONTIG_HASHES_MATCH } from '../modules/local/verify_contig_hashes_match/main'
include { FORMAT_OUTPUT_RESULTS      } from '../modules/local/format_output_results/main'

workflow MAG_ASSEMBLY_LINKING_PIPELINE {
    main:
        processed_acc_ch = params.processed_acc ? Channel.fromPath(params.processed_acc, checkIfExists: true) : []

        // If custom input accessions are provided, use them instead of the downloaded accessions
        if (params.external_input) {
            accessions_list_ch = Channel.fromPath(params.external_input, checkIfExists: true)
            metadata_ch = params.external_metadata ? Channel.fromPath(params.external_metadata, checkIfExists: true) : []

        } else {
            COLLECT_INPUT_ACCESSIONS(processed_acc_ch, params.input_accessions, params.gut_mapping, params.catalogue_metadata)
            metadata_ch = COLLECT_INPUT_ACCESSIONS.output.metadata
            accessions_list_ch = COLLECT_INPUT_ACCESSIONS.output.input_accessions
        }

        // Input accessions are splitted to process them faster in parallel tasks
        accessions_batches_ch = accessions_list_ch.splitText(by: params.batch_size, file: "batch")

        MAP_GENOMES_TO_ASSEMBLIES(accessions_batches_ch, metadata_ch, params.catalogue_metadata, params.gut_mapping)

        VERIFY_CONTIG_HASHES_MATCH(MAP_GENOMES_TO_ASSEMBLIES.OUT)

        mag_assembly_pairs_ch = VERIFY_CONTIG_HASHES_MATCH.output.mag_assembly_pairs
        not_linked_mags_ch = VERIFY_CONTIG_HASHES_MATCH.output.not_linked_mags
        previous_table_ch = params.previous_table ? Channel.fromPath(params.previous_table) : []

        FORMAT_OUTPUT_RESULTS(mag_assembly_pairs_ch.collect(), not_linked_mags_ch.collect(), metadata_ch, processed_acc_ch, previous_table_ch)

    emit:
        mag_to_assembly_links_ch = FORMAT_OUTPUT_RESULTS.out.mag_to_assembly_links
        processed_accessions_ch = FORMAT_OUTPUT_RESULTS.out.processed_accessions
}
