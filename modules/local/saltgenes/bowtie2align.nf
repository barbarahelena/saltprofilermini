process SALTGENES_BOWTIE2ALIGN {
    tag "$subjectid"
    label 'process_medium'

    conda "bioconda::bowtie2=2.4.2 bioconda::samtools=1.11 conda-forge::pigz=2.3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0':
        'biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0' }"

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
    INDEX=\$(find -L ./ -name "*.rev.1.bt2l" -o -name "*.rev.1.bt2" | sed 's/.rev.1.bt2l//' | sed 's/.rev.1.bt2//')
    bowtie2 \\
        -p "${task.cpus}" \\
        -x \$INDEX \\
        $args \\
        -1 "${reads[0]}" -2 "${reads[1]}" \\
        1> /dev/null \\
        2> ${subjectid}.bowtie2.log > bowtieout.sam
    
    ls -lh bowtieout.sam

    samtools view -@ "${task.cpus}" -bS bowtieout.sam > output.bam
    samtools sort -@ "${task.cpus}" -o ${subjectid}.bam output.bam
    samtools index ${subjectid}.bam

    samtools idxstats ${subjectid}.bam > counts_${subjectid}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${subjectid}"

    """
    touch ${subjectid}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
