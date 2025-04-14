process SALTGENES_MERGE {
    label "process_single"

    conda "bioconda::csvkit=1.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/csvkit:2.1.0--2074e44ef35da83c' :
        'community.wave.seqera.io/library/csvkit:2.1.0--89933f91a588cea2' }"

    input:
    path(summaries)

    output:
    path("*.csv")      , emit: combined
    path "versions.yml", emit: versions

    script:
    """
    csvstack ${summaries} > all_saltgenes.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvkit: \$(csvstack --version | head -n 1)
    END_VERSIONS
    """
}
