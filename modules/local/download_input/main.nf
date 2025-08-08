process DOWNLOAD_INPUT_BINS_AND_MAGS {
    container 'quay.io/microbiome-informatics/mag-assembly-linking:v1.1'
    label 'process_single'

    input:
    path processed_acc
    val input_accessions
    path gut_mapping
    val catalogue_metadata


    output:
    path "${input_accessions}", emit: input_accessions
    path catalogue_metadata, emit: metadata

    script:
    def processed_accession = processed_acc ? "--processed-acc ${processed_acc}" : ""
    """
    download_genome_accessions.py \\
        $processed_accession \\
        --output-accessions ${input_accessions} \\
        --gut-mapping ${gut_mapping} \\
        --catalogue-metadata ${catalogue_metadata}
    """
}
