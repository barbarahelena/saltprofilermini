process SALTGENES_CATPERGENE {
    tag "$subjectid, $gene"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(subjectid), val(binid), val(gene), path(fasta)
    tuple val(subjectid), val(binid), val(gene), path(gff)

    output:
    tuple val(subjectid), val(gene), path("*.fasta")  , emit: mergedfa
    tuple val(subjectid), val(gene), path("*.gff")    , emit: mergedgff
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}"

    """
    cat $fasta  > ${subjectid}_${gene}.fasta
    cat $gff  > ${subjectid}_${gene}.gff

    cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}"

    """
    touch ${subjectid}_${gene}.fasta

    cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
