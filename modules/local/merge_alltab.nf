process MERGE_ALLTAB {

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    path(saltgenes)
    path(tax)

    output:
    path("saltgenes_tax_depths.tsv"  )   , emit: summary
    path "versions.yml"                  , emit: versions

    script:
    """
    merge_alltab.py --saltgenes ${saltgenes} --tax ${tax} --out saltgenes_tax_depths.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
