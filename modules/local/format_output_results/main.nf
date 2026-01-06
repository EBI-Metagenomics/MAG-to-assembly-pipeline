process FORMAT_OUTPUT_RESULTS {
    label 'process_single'

    container 'quay.io/microbiome-informatics/mag-assembly-linking:v1.1-multiarch'

    input:
    path mapped_genomes_tsv
    path not_mapped_genomes_tsv
    path catalogues_metadata
    path previously_processed_acc, stageAs: 'previously_processed_acc.tsv'
    path previous_table, stageAs: 'previous_table.tsv'

    output:
    path 'mag_to_assembly_mapping_*.tsv'       , emit: mag_to_assembly_mapping
    path 'processed_accessions_*.tsv'          , emit: processed_accessions
    path 'mag_to_assembly_links_to_unlink*.tsv', emit: links_to_unlink, optional: true

    script:
    def processed_accessions   = previously_processed_acc ? "previously_processed_acc.tsv" : ""
    // If proccessing was done from scratch and previous table is provided, a file will be created
    // that lists MAG-to-assembly links that were in the previous table but not in the new results
    def write_deleted_flag     = (!previously_processed_acc && previous_table) ? "--write-deleted" : ""
    def previous_linking_table = previous_table ? "--previous-table previous_table.tsv" : ""
    def metadata_file          = catalogues_metadata ? "--catalogue-metadata ${catalogues_metadata}" : ""

    """
    # Combine all processed accessions (mapped + not mapped + previously processed)
    cat ${mapped_genomes_tsv} ${not_mapped_genomes_tsv} ${processed_accessions} | cut -f 1 > processed_accessions_\$(date +"%Y-%m-%d_%Hh%Mm").tsv

    format_output.py \\
        ${previous_linking_table} \\
        ${write_deleted_flag} \\
        ${metadata_file}  \\
        ${mapped_genomes_tsv}
    """

    stub:
    """
    touch mag_to_assembly_mapping_\$(date +"%Y-%m-%d_%Hh%Mm").tsv
    touch processed_accessions_\$(date +"%Y-%m-%d_%Hh%Mm").tsv
    touch mag_to_assembly_links_to_unlink\$(date +"%Y-%m-%d_%Hh%Mm").tsv
    """
}
