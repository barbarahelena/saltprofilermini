process SALTGENES_CATPERASSEMBLY {
    tag "$meta.id"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), val(gene), path(fasta), path(gff)

    output:
    tuple val(meta), val(gene), path("*_allgenes.csv")  , emit: merged
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def output_csv = "${meta.id}_allgenes.csv"

    """
    # Write the header to the output CSV
    echo "Assembly,Gene,ProteinID,Sequence" > "${output_csv}"

    # Loop through all GFF and FASTA files
    for gff_file in *.gff; do
        # Find files
        base_name=\$(basename "\$gff_file" .gff)
        fasta_file="\${base_name}_fixed.fasta"

        # Check if the GFF file contains valid data
        if ! grep -q "gene=" "\$gff_file"; then
            echo "Skipping \$gff_file: No valid gene data found."
            continue
        fi

        # Extract the assembly name
        assembly=\$(echo "\$base_name" | cut -d'_' -f1-2)

        # Extract the assembly name, gene name, and set number from the GFF file
        while IFS=\$'\\t' read -r seqid start end strand attributes; do
            # Skip lines that are not CDS features
            if [[ "\$attributes" != *"gene="* ]]; then
                continue
            fi

            # Parse the attributes field to extract gene and other information
            gene=\$(echo "\$attributes" | grep -oP "gene=\\K[^;]+")
            other_data=\$(echo "\$attributes" | grep -oP "ID=\\K[^;]+")

            # Extract the sequence from the FASTA file
            sequence=\$(grep -A1 "\$seqid" "\$fasta_file" | tail -n1)

            echo "Processing gene: \$gene"

            # Skip if the sequence is empty
            if [[ -z "\$sequence" ]]; then
                echo "Skipping \$gff_file: No sequence found for \$gene."
                continue
            fi

            # Write the data to the output CSV
            echo "\$assembly,\$gene,\$other_data,\$sequence" >> "${output_csv}"
        done < "\$gff_file"
    done

    echo "Merged results saved to ${output_csv}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version | sed 's/^.*coreutils) //; s/ .*\$//')
        bash: \$(bash --version | head -n 1 | awk '{print \$4}')
        grep: \$(grep --version | head -n 1 | awk '{print \$4}')
    END_VERSIONS
    """

    stub:
    def output_csv = "${meta.id}_allgenes.csv"

    """
    touch ${output_csv}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version | sed 's/^.*coreutils) //; s/ .*\$//')
        bash: \$(bash --version | head -n 1 | awk '{print \$4}')
        grep: \$(grep --version | head -n 1 | awk '{print \$4}')
    END_VERSIONS
    """
}