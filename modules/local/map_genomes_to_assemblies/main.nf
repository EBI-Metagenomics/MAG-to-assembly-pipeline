process MAP_GENOMES_TO_ASSEMBLIES {
    container 'quay.io/microbiome-informatics/mag-assembly-linking:v1.1'

    input:
    path accessions

    output:
    path "*.tsv", emit: genome_to_assemblies_mapping

    script:
    def debug_flag = params.debug ? "--debug" : ""
    """
    map_genomes_to_assemblies.py \\
        -i ${accessions} \\
        -o ${accessions}.links.tsv \\
        -p ${accessions}.not_linked.tsv \\
        -f ${accessions}.failed.tsv \\
        ${debug_flag}
    """

    stub:
    """

    """
}
