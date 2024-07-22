process FIND_PRIMARY_ASSEMBLY {
    input:
    path input

    output:
    path "*.confirmed.tsv"

    script:
    """
    link_MAG_to_primary_metagenome_assembly.py -i ${input} -o ${input}.confirmed.tsv -p ${input}.putative.tsv -f ${input}.failed.tsv --download-folder fastas
    """
}