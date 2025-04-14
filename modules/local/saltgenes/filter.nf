process SALTGENES_FILTER {
    tag "$meta.id"
    label "process_single"

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2':
        'biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"

    input:
    tuple val(meta), path(fasta), path(gff), val(gene)

    output:
    tuple val(meta), val(gene), path("*_${gene}_fixed.fasta"), path("*_${gene}.gff")   , emit: seqs
    path "versions.yml"                                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "Processing GFF: ${gff} for gene: ${gene}"

    # Extract relevant entries from the GFF file
    awk -F '\\t' -v gene="${gene}" '
    BEGIN { OFS="\\t"; IGNORECASE=1 }
    \$3 == "CDS" && tolower(\$9) ~ "gene="tolower(gene)"([;]|\$)" {
        print \$1, \$4, \$5, \$7, \$9
    }
    ' "${gff}" > "${meta.id}_${gene}.gff"

    echo "Extracted gene coordinates to ${meta.id}_${gene}.gff"

    # Check if the gene was found
    if [ -s "${meta.id}_${gene}.gff" ]; then
        echo "Gene ${gene} found in ${gff}"
    else
        echo "Gene ${gene} not found in ${gff}"
    fi

    # Extract the sequence using bedtools
    if [ -s "${meta.id}_${gene}.gff" ]; then
        bedtools getfasta \\
            -fi "${fasta}" \\
            -bed "${meta.id}_${gene}.gff" \\
            -s \\
            -fo "${meta.id}_${gene}_fixed.fasta"
        echo "Extracted sequences to ${meta.id}_${gene}_fixed.fasta"
    else
        echo "No sequences extracted for ${gene}."
        touch "${meta.id}_${gene}_fixed.fasta"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version 2>&1 | head -n 1 | awk '{print \$3}')
        bedtools: \$(bedtools --version | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${gene}_fixed.fasta
    touch ${meta.id}_${gene}.gff3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$3}')
        bedtools: \$(bedtools --version | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """
}
