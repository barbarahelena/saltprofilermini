process SALTGENES_BWAMEM2BUILD {
    tag "$subjectid"

    conda "bioconda::bwa-mem2=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa-mem2%3A2.2.1--he513fc3_0':
        'biocontainers/bwa-mem2:2.2.1--he513fc3_0' }"

    input:
    tuple val(subjectid), path(fasta), path(gff)

    output:
    tuple val(subjectid), path("${subjectid}_bwamem2"), path(gff)  , emit: index
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}"

    """
    mkdir "${subjectid}_bwamem2"
    bwa-mem2 \\
        index \\
        $fasta \\
        $args \\
        -p ${subjectid}_bwamem2/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(echo \$(bwa-mem2 --version 2>&1) | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}"
    """
    mkdir "${subjectid}_bwamem2"
    touch ${subjectid}_bwamem2/${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(echo \$(bwa-mem2 --version 2>&1) | sed 's/^.* //')
    END_VERSIONS
    """
}
