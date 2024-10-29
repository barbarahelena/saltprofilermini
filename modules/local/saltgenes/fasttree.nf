process SALTGENES_FASTTREE {
    tag "${subjectid}, ${gene}"
    label 'process_single'

    conda "bioconda::fasttree=2.1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fasttree:2.1.11--h031d066_4':
        'biocontainers/fasttree:2.1.11--h031d066_4' }"

    input:
    tuple val(subjectid), val(gene), path(msa)

    output:
    tuple val(subjectid), val(gene), path("*.tree")   , emit: tree
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}"

    """
    fasttree \\
        -nt \\
        -gtr \\
        -gamma \\
        $msa \\
        > ${subjectid}_${gene}.tree

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: \$(fasttree -help 2>&1 | head -n 1 | sed -n 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}_${gene}"
    """
    touch ${subjectid}_${gene}.tree

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: \$(fasttree -help 2>&1 | head -n 1 | sed -n 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
