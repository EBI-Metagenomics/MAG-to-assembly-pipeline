process POSTPROCESSING {
    input:
    path input
    path catalogue_metadata

    output:
    path 'RETROFIT.tsv'

    script:
    """
    finalise_output.py --previous-table ${params.previous_table} --catalogue-metadata ${catalogue_metadata} ${input}
    """
}