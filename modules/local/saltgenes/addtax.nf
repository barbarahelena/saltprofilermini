process SALTGENES_ADDTAX {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bioconductor-ape bioconda::r-dplyr bioconda::ggplot2 bioconda::r-tidyr bioconda::r-forcats bioconda::r-ggthemes"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(bam), path(tax)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    

    writeLines(c("\\"${task.process}\\":", 
        paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
        paste0("    ggtree: ", packageVersion("ggtree")) ),
        "versions.yml")
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        saltgenes: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
