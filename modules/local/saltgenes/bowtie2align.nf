process SALTGENES_BOWTIE2ALIGN {
    tag "$subjectid, $gene"
    label 'process_medium'

    conda "bioconda::bowtie2=2.4.2 bioconda::samtools=1.11 conda-forge::pigz=2.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0':
        'biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0' }"

    input:
    tuple val(subjectid), val(gene), path(index), path(reads), path(gff)

    output:
    tuple val(subjectid), val(gene), path("${prefix}.bam"), path(gff)  , emit: bam
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}_${gene}"
    def forwardreads = reads[0]
    def revreads = reads[1]

    """
    bowtie2 \\
        -p "${task.cpus}" \\
        -x ${prefix} \\
        -S /dev/null \\
        $args \\
        -1 ${forwardreads} \\
        -2 ${revreads} \\
        2> "${prefix}.bowtie2.log" | \
        samtools view -@ "${task.cpus}" -bS | \
        samtools sort -@ "${task.cpus}" -o "${prefix}.bam"
    samtools index "${prefix}.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}_${gene}"

    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
