process SALTGENES_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2':
        'biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"

    input:
    tuple val(meta), path(gff), path(fasta), val(gene)

    output:
    tuple val(meta), val(gene), path("*.fasta")     , emit: seqs
    tuple val(meta), val(gene), path("*.gff")       , emit: gff
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${gene}"

    """
    awk '\$3 == "gene" && \$0 ~ /Name=${gene}/' ${gff} > ${prefix}.gff

    bedtools \\
        getfasta \\
        $args \\
        -bed ${prefix}.gff \\
        -fi ${fasta} \\
        -fo ${prefix}.fasta
    
    awk -v sample_id="bin_${meta.id}" -v gene_id="gene_${gene}" '/^>/ {print ">" sample_id "_" gene_id "_" substr($0, 2)} !/^>/ {print}' \\
            ${prefix}.fasta > ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version | sed -n 's/.*version \\([0-9.]*\\)/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}/${prefix}_${gene}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread: \$(gffread --version | sed -n 's/.*version \\([0-9.]*\\)/\\1/p')
    END_VERSIONS
    """
}
