process DOWNLOAD_INPUT {
    label "python_based"
    conda "${moduleDir}/environment.yml"
    
    input:
    path processed_acc
    val input_accessions
    path gut_mapping
    val catalogue_metadata
    

    output:
    path "${input_accessions}", emit: input_accessions
    path catalogue_metadata, emit: metadata

    script:
    """
    download_genome_accessions.py --processed-acc ${processed_acc} --output-accessions ${input_accessions} --gut-mapping ${gut_mapping} --catalogue-metadata ${catalogue_metadata}
    """
}