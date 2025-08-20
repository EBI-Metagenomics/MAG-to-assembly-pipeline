process VERIFY_CONTIG_HASHES_MATCH {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/microbiome-informatics/mag-assembly-linking:v1.1'

    input:
    tuple val(meta), path(genome_to_assemblies_mapping)

    output:
    tuple val(meta), path("*.verified.tsv"), emit: verified_pairs, optional: true
    path("*.invalid.tsv")                  , emit: invalid_pairs, optional: true
    path("*.err")                          , emit: error_log, optional: true

    script:
    def prefix = task.ext.prefix ?: "$meta.id"
    def cleanup_flag = params.cleanup ? "--cleanup" : ""
    def debug_flag = params.debug ? "--debug" : ""

    """
    verify_contig_hashes_match.py \\
        --input ${genome_to_assemblies_mapping} \\
        --output_verified ${prefix}.verified.tsv \\
        --output_invalid ${prefix}.invalid.tsv \\
        --errors ${prefix}.err \\
        ${cleanup_flag} \\
        ${debug_flag}
    """

    stub:
    def prefix = task.ext.prefix ?: "$meta.id"

    """
    touch "${prefix}.verified.tsv"
    touch "${prefix}.invalid.tsv"
    touch "${prefix}.err"
    """
}
