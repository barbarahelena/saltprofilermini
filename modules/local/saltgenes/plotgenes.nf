process SALTGENES_PLOTGENES {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bioconductor-ape bioconda::r-dplyr bioconda::ggplot2 bioconda::r-tidyr bioconda::r-forcats bioconda::r-ggthemes"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(bam)

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
    tre <- ape::read_tree("HELIBA_100018_galE.tree")
    gff <- read.delim("HELIBA_100018_galE.gff", sep = "\t", header = FALSE)
    # taxonomy ? 

    bin <- str_extract(tre\$tip.label, "HELI[A-Z]*_[0-9]*.[0-9]*")
    metadata <- data.frame(node = 1:length(tre\$tip.label), 
                            bin = bin,
                            label = tre\$tip.label)

    pl <- ggtree(tre) %<+% metadata + geom_label(aes(label = bin)) + ggtitle("HELIBA_100018, galE")
    ggsave("testtree.pdf", pl, width = 16, height = 6)

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
