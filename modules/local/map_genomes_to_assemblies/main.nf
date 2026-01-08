process MAP_GENOMES_TO_ASSEMBLIES {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/microbiome-informatics/mag-assembly-linking:v1.1-multiarch'

    input:
    tuple val(meta), path(accessions)

    output:
    tuple val(meta), path("*.mapped_pairs.tsv"), emit: tsv_mapping, optional: true
    path("*.no_assembly.tsv")                  , emit: no_assembly_found, optional: true
    path("*.err")                              , emit: error_log, optional: true

    script:
    def prefix     = task.ext.prefix ?: "$meta.id"
    def debug_flag = params.debug ? "--debug" : ""

    """
    map_genomes_to_assemblies.py \\
        --input ${accessions} \\
        --output ${prefix}.mapped_pairs.tsv \\
        --no_assembly_found ${prefix}.no_assembly.tsv \\
        --errors ${prefix}.err \\
        ${debug_flag}
    """

    stub:
    def prefix = task.ext.prefix ?: "$meta.id"

    """
    touch "${prefix}.mapped_pairs.tsv"
    touch "${prefix}.no_assembly.tsv"
    touch "${prefix}.err"
    """
}
