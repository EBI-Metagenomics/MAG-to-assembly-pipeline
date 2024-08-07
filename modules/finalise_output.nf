process POSTPROCESSING {
    label "python_based"
    conda "${moduleDir}/environment.yml"
    
    input:
    path input
    path not_linked_mags
    path catalogue_metadata
    // inputs below are optional and can be placeholder empty files 
    // it is required to prevent nexflow from crashing because by default it does not support optional inputs
    path previous_processed_acc, stageAs: 'previous_processed_acc.tsv'
    path previous_table, stageAs: 'previous_table.tsv'


    output:
    path 'mag_to_assembly_links_*.tsv'
    path 'processed_accessions_*.tsv'

    script:
    """
    
    cat ${input} ${not_linked_mags} previous_processed_acc.tsv | cut -f 1 > processed_accessions_\$(date +"%Y-%m-%d_%Hh%Mm").tsv

    finalise_output.py \
        --previous-table previous_table.tsv \
        --catalogue-metadata ${catalogue_metadata} \
        ${input}

    """
}