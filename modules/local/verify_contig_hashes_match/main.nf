process VERIFY_CONTIG_HASHES_MATCH {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/microbiome-informatics/mag-assembly-linking:v1.1'

    input:
    tuple val(meta), path(genome_to_assemblies_mapping)

    output:
    tuple val(meta), path("*.validated.tsv"), emit: validated_mag_assembly_pairs
    tuple val(meta), path("*.err")          , emit: error_log

    script:
    def prefix = task.ext.prefix ?: "$meta.id"
    def cleanup_flag = params.cleanup ? "--cleanup" : ""
    def debug_flag = params.debug ? "--debug" : ""

    """
    verify_contig_hashes_match.py \\
        --input ${genome_to_assemblies_mapping} \\
        --output ${prefix}.validated.tsv \\
        --errors ${prefix}.err \\
        ${cleanup_flag} \\
        ${debug_flag}
    """

    stub:
    def prefix = task.ext.prefix ?: "$meta.id"

    """
    touch "${prefix}.validated.tsv"
    touch "${prefix}.err"
    """
}
