process SALTGENES_BOWTIE2BUILD {
    tag "$subjectid"
    label 'process_medium'

    conda "bioconda::bowtie2=2.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.5.4--h7071971_4':
        'biocontainers/bowtie2:2.5.4--h7071971_4' }"

    input:
    tuple val(subjectid), path(fasta), path(gff)

    output:
    tuple val(subjectid), path("*.bt2"), path(gff)  , emit: index
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}"

    """
    mkdir bowtie
    bowtie2-build --threads $task.cpus $fasta $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}"
    """
    touch ${prefix}.bt2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
