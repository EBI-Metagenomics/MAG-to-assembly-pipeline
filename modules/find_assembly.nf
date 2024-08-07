process FIND_PRIMARY_ASSEMBLY {
    label "python_based"
    conda "${moduleDir}/environment.yml"
    
    input:
    path input

    output:
    path "*.links.tsv", emit: mag_assembly_pairs
    path "*.not_linked.tsv", emit: not_linked_mags

    script:
    """
    
    link_MAG_to_primary_metagenome_assembly.py \
        -i ${input} \
        -o ${input}.links.tsv \
        -p ${input}.putative.not_linked.tsv \
        -f ${input}.failed.not_linked.tsv \
        --download-folder fastas \
        --cleanup

    """
}