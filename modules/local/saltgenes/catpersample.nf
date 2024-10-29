process SALTGENES_CATPERSAMPLE {
    tag "${subjectid}"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(subjectid), val(gene), path(fasta), path(gff)

    output:
    tuple val(subjectid), path("*.fasta"), path("*.gff")  , emit: mergedgenes
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}"

    """
    cat $fasta  > ${subjectid}.fasta
    cat $gff >> ${subjectid}.gff

    cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}"

    """
    touch ${subjectid}.fasta
    touch ${subjectid}.gff

    cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
