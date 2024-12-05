process FIND_PRIMARY_ASSEMBLY {
    container 'quay.io/microbiome-informatics/mag-assembly-linking:v1.1'
    
    input:
    path accessions

    output:
    path "*.links.tsv", emit: mag_assembly_pairs
    path "*.not_linked.tsv", emit: not_linked_mags

    script:
    def cleanup_flag = params.cleanup ? "--cleanup" : ""
    def debug_flag = params.debug ? "--debug" : ""
    """
    link_MAG_to_primary_metagenome_assembly.py \\
        -i ${accessions} \\
        -o ${accessions}.links.tsv \\
        -p ${accessions}.putative.not_linked.tsv \\
        -f ${accessions}.failed.not_linked.tsv \\
        --download-folder fastas \\
        ${cleanup_flag} \\
        ${debug_flag} 
    """
}