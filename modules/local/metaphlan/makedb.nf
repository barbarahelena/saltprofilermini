process METAPHLAN_MAKEDB {
    label 'process_medium'
    storeDir 'db'

    conda "bioconda::metaphlan=4.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaphlan:4.0.5--pyhca03a8a_0' :
        'biocontainers/metaphlan:4.0.5--pyhca03a8a_0' }"

    output:
    path "metaphlan_db"         , emit: db

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    metaphlan \\
        --install \\
        --nproc $task.cpus \\
        --bowtie2db metaphlan_db \\
        $args
    """
}
