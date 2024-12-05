process MERGE_TAXDEPTHS {

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    path(depths)
    path(tax)

    output:
    path("bins_taxdepths.tsv")  , emit: summary
    path "versions.yml"         , emit: versions

    script:
    """
    mergetax_depths.py --depths_${depths} --tax ${tax} --out bin_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
