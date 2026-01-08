process COLLECT_INPUT_ACCESSIONS {
    label 'process_single'

    container 'quay.io/microbiome-informatics/mag-assembly-linking:v1.1-multiarch'

    input:
    path skip_accessions
    path gut_accessions_mapping

    output:
    path "input_accessions.tsv"    , emit: input_accessions
    path "all_catalog_metadata.tsv", emit: catalogues_metadata

    script:
    def skip_accessions_flag = skip_accessions ? "--skip_accessions ${skip_accessions}" : ""

    """
    download_genome_accessions.py \\
        $skip_accessions_flag \\
        --gut-mapping ${gut_accessions_mapping} \\
        --output-accessions input_accessions.tsv \\
        --catalogue-metadata all_catalog_metadata.tsv
    """

    stub:
    """
    touch input_accessions.tsv
    touch all_catalog_metadata.tsv
    """
}
