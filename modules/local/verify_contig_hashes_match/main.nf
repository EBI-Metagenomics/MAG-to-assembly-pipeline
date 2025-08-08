process VERIFY_CONTIG_HASHES_MATCH {
    container 'quay.io/microbiome-informatics/mag-assembly-linking:v1.1'

    input:
    path genome_to_assemblies_mapping

    output:
    path "*.tsv", emit: validated_mag_assembly_pairs
    path "*.err", emit: error_log

    script:
    def cleanup_flag = params.cleanup ? "--cleanup" : ""
    def debug_flag = params.debug ? "--debug" : ""
    """
    verify_contig_hashes_match.py \\
        -i ${genome_to_assemblies_mapping} \\
        -o ${genome_to_assemblies_mapping}.validated.tsv \\
        --download-folder fastas \\
        ${cleanup_flag} \\
        ${debug_flag}
    """

    stub:
    """

    """
}
