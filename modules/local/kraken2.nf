process KRAKEN2 {
    tag "${sample}"

    conda "bioconda::kraken2=2.0.8_beta"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kraken2:2.0.8_beta--pl526hc9558a2_2' :
        'biocontainers/kraken2:2.0.8_beta--pl526hc9558a2_2' }"

    input:
    tuple val(sample), path(bins)
    path(database)
    path(taxonomy)

    output:
    tuple val(sample), path("${sample}_taxonomy.tsv")     , emit: tax
    tuple val(sample), path("${sample}_headers.txt")      , emit: headers
    path "versions.yml"                                   , emit: versions

    script:
    prefix = task.ext.prefix ?: "${sample}"

    """
    set -e  # Exit on error
    set -u  # Treat unset variables as an error
    set -o pipefail  # Exit if any command in a pipeline fails

    output_file="${sample}_taxonomy.tsv"
    echo -e "Sample\tBin\tTaxID\tTaxonomy\tContigsMatch\tContigsTotal\tComments" > "\$output_file"

    header_file="${sample}_headers.txt"

    echo "Processing bins: ${bins}"
    for bin in ${bins}; do
        echo "Processing bin: \$bin"
        binname=\$(basename "\$bin" .fa)
        echo "Bin name: \$binname"
        
        echo "Running Kraken2 on \$binname"
        kraken2 --db ${database} \\
                --threads ${task.cpus} \\
                "\$bin" \\
                --use-names \\
                --report "kraken2_\${binname}_report.txt" \\
                --classified-out "out_\${binname}_classified.fa"

        echo "Done with Kraken2, now making the table.."
        cat "out_\${binname}_classified.fa" | grep '>' | cut -d'|' -f 2 | sort | uniq -c | sort -nrk1 | head -n10
        cat "out_\${binname}_classified.fa" | grep '>' >> \$header_file

        sampleid="${sample}"
        taxid_info=\$(cat "out_\${binname}_classified.fa" | grep '>' | cut -d'|' -f 2 | sort | uniq -c | sort -nrk1 | head -n2)
        first_taxid=\$(echo "\$taxid_info" | head -n1 | awk '{print \$2}')
        second_taxid=\$(echo "\$taxid_info" | tail -n1 | awk '{print \$2}')
        
        if [ "\$first_taxid" = "2" ]; then
            taxid="\$second_taxid"
            contigmatch=\$(cat "out_\${binname}_classified.fa" | grep '>' | cut -d'|' -f 2 | sort | uniq -c | sort -nrk1 | head -n2 | tail -n1 | awk '{print \$1}')
            bacmatch=\$(cat "out_\${binname}_classified.fa" | grep '>' | cut -d'|' -f 2 | sort | uniq -c | sort -nrk1 | head -n1 | awk '{print \$1}')
            comment="Most common classification was Bacteria with \$bacmatch matches"
            echo "First TaxID was 2, using second most common TaxID: \$taxid"
        else
            taxid="\$first_taxid"
            contigmatch=\$(cat "out_\${binname}_classified.fa" | grep '>' | cut -d'|' -f 2 | sort | uniq -c | sort -nrk1 | head -n1 | awk '{print \$1}')
            comment=""
            echo "Using most common TaxID: \$taxid"
        fi
        
        echo "Final taxID: \$taxid"
        
        echo "Extracting taxonomy"
        if [ -z "\$taxid" ]; then
            echo "Warning: TaxID is empty. Setting taxonomy to 'Unknown'"
            taxonomy_name="Unknown"
        else
            # Use awk to find exact match and extract the second column
            taxonomy_name=\$(awk -F'\t' -v tid="\$taxid" '\$1 == tid {for (i = 2; i <= NF; i++) printf "%s ", \$i; print ""; exit}' "${taxonomy}")
            if [ -z "\$taxonomy_name" ]; then
                echo "Warning: No exact match found for TaxID \$taxid. Setting taxonomy to 'Unknown'"
                taxonomy_name="Unknown"
            fi
        fi
        echo "Taxonomy: \$taxonomy_name"

        echo "Calculating contigs that match this taxonomy, and total number of contigs"
        
        contigtotal=\$(cat "out_\${binname}_classified.fa" | grep '>' | wc -l)
        echo "Contigs match: \$contigmatch, Total contigs: \$contigtotal"

        echo "Appending to output file"
        echo -e "\${sampleid}\t\${binname}\t\${taxid}\t\${taxonomy_name}\t\${contigmatch}\t\${contigtotal}\t\${comment}" >> "\$output_file"
    done
    
    echo "Contents of output file:"
    cat "\$output_file"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //' | sed 's/ Copyright.*//')
    END_VERSIONS
    """
}
