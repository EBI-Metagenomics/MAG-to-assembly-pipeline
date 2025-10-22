#!/usr/bin/env nextflow

include { COLLECT_INPUT_ACCESSIONS   } from '../modules/local/collect_input_accessions/main'
include { MAP_GENOMES_TO_ASSEMBLIES  } from '../modules/local/map_genomes_to_assemblies/main'
include { VERIFY_CONTIG_HASHES_MATCH } from '../modules/local/verify_contig_hashes_match/main'
include { FORMAT_OUTPUT_RESULTS      } from '../modules/local/format_output_results/main'

workflow MAG_ASSEMBLY_LINKING_PIPELINE {
    main:
        skip_accessions_ch = params.skip_accessions ? Channel.fromPath(params.skip_accessions, checkIfExists: true) : []

        // If custom list of input accessions is provided, use it instead of the accessions collected from ENA and MGnify
        if (params.accessions_list) {
            accessions_list_ch     = Channel.fromPath(params.accessions_list, checkIfExists: true)
            catalogues_metadata_ch = params.catalogues_metadata ? Channel.fromPath(params.catalogues_metadata, checkIfExists: true) : []

        // Otherwise, build list of input genomes from ENA bins and MAGs and MGnify catalogues
        } else {
            COLLECT_INPUT_ACCESSIONS(skip_accessions_ch, params.gut_accessions_mapping)
            accessions_list_ch     = COLLECT_INPUT_ACCESSIONS.output.input_accessions
            catalogues_metadata_ch = COLLECT_INPUT_ACCESSIONS.output.catalogues_metadata
        }

        // Input accessions are splitted into batches of size params.batch_size to process them faster in parallel
        accessions_batches_ch = accessions_list_ch
            .splitText(by: params.batch_size, file: "batch")
            .map { batch_file ->
                def meta = [id: batch_file.name]
                [meta, batch_file]
            }

        // Find primary assembly for each genome using information from ENA
        MAP_GENOMES_TO_ASSEMBLIES(accessions_batches_ch)

        // Verify that contigs are identical in a genome and its assembly
        VERIFY_CONTIG_HASHES_MATCH(MAP_GENOMES_TO_ASSEMBLIES.output.tsv_mapping)

        mag_assembly_pairs_ch = VERIFY_CONTIG_HASHES_MATCH.output.verified_pairs
            .map { _meta, pairs_tsv -> pairs_tsv}
            .collectFile(name: "mag_to_assembly_mapping.tsv")
        no_assembly_genomes_ch = VERIFY_CONTIG_HASHES_MATCH.output.invalid_pairs
            .mix(MAP_GENOMES_TO_ASSEMBLIES.output.no_assembly_found)
            .collectFile(name: "no_assembly_genomes.tsv")
        previous_results_ch = params.merge_with_results ? Channel.fromPath(params.merge_with_results) : []

        // Format output results: create a table with MAGs, their primary assemblies and MGYG accessions,
        // and update the list of processed accessions
        FORMAT_OUTPUT_RESULTS(
            mag_assembly_pairs_ch,
            no_assembly_genomes_ch,
            catalogues_metadata_ch,
            skip_accessions_ch,
            previous_results_ch
        )

    emit:
        mag_to_assembly_mapping = FORMAT_OUTPUT_RESULTS.out.mag_to_assembly_mapping
        processed_accessions    = FORMAT_OUTPUT_RESULTS.out.processed_accessions
}
