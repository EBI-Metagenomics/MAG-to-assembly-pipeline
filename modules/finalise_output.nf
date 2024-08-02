process POSTPROCESSING {
    label "python_based"
    conda "${moduleDir}/environment.yml"
    
    input:
    path input
    path catalogue_metadata
    path previous_table

    output:
    path 'mag_to_assembly_links.tsv'

    script:
    """
    finalise_output.py --previous-table ${previous_table} --catalogue-metadata ${catalogue_metadata} ${input}
    """
}