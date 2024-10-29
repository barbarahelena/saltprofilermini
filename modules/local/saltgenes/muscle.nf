process SALTGENES_MUSCLE {
    tag "$subjectid, $gene"
    label 'process_single'

    conda "bioconda::mafft=7.525"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/muscle:3.8.1551--h2d50403_3':
        'biocontainers/muscle:3.8.1551--h2d50403_3' }"

    input:
    tuple val(subjectid), val(gene), path(fasta), path(gff)

    output:
    tuple val(subjectid), val(gene), path("*_aln.fasta")  , emit: msa
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    muscle -in ${fasta} -out ${subjectid}_${gene}_aln.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS

    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch ${subjectid}_${gene}_aln.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}