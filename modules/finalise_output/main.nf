process FINALISE_OUTPUT {
    label "python_based"
    conda "${moduleDir}/environment.yml"
    
    input:
    path linked_tsv
    path not_linked_tsv
    path catalogue_metadata
    path previous_processed_acc, stageAs: 'previous_processed_acc.tsv'
    path previous_table, stageAs: 'previous_table.tsv'


    output:
    path 'mag_to_assembly_links_*.tsv'
    path 'processed_accessions_*.tsv'

    script:
    def processed_accessions = previous_processed_acc ? "previous_processed_acc.tsv" : ""
    def previous_linking_table = previous_table ? "--previous-table previous_table.tsv" : ""
    def metadata_file = catalogue_metadata ? "--catalogue-metadata ${catalogue_metadata}" : ""
    """

    cat ${linked_tsv} ${not_linked_tsv} ${processed_accessions} | cut -f 1 > processed_accessions_\$(date +"%Y-%m-%d_%Hh%Mm").tsv

    finalise_output.py \
        ${previous_linking_table} \
        ${metadata_file}  \
        ${input}

    """
}