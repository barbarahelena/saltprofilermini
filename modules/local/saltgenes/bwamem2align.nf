process SALTGENES_BWAMEM2ALIGN {
    tag "$subjectid"
    label 'process_medium'

    conda "bioconda::bwa-mem2=2.2.1,bioconda::htslib=1.19.1,bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0':
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0' }"

    input:
    tuple val(subjectid), path(index), path(gff), path(reads)

    output:
    tuple val(subjectid), path("${subjectid}.bam"), path(gff)  , emit: bam
    tuple val(subjectid), path("*.txt")                        , emit: stats
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}"

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    bwa-mem2 mem \\
        -t "${task.cpus}" \\
        \$INDEX \\
        $args \\
        "${reads[0]}" "${reads[1]}" \\
        1> /dev/null \\
        2> ${subjectid}.bwa-mem2.log > bwaout.sam

    ls -lh bwaout.sam

    samtools view -@ "${task.cpus}" -bS bwaout.sam -F 256 -q 5 > output.bam
    samtools sort -@ "${task.cpus}" -o ${subjectid}.bam output.bam
    samtools index ${subjectid}.bam

    samtools idxstats ${subjectid}.bam > counts_${subjectid}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(echo \$(bwa-mem2 --version 2>&1) | sed 's/^.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}"

    """
    touch ${subjectid}.bam
    touch counts_${subjectid}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(echo \$(bwa-mem2 --version 2>&1) | sed 's/^.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
