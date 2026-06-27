process MOSDEPTH_HIGH_COVERAGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gawk:5.3.0'
        : 'biocontainers/gawk:5.3.0'}"

    input:
    tuple val(meta), path(regions_bed)

    output:
    tuple val(meta), path("*.high_coverage.bed"), emit: bed
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def threshold = task.ext.threshold ?: 50000
    // mosdepth regions BED (--by 500): col4 is the mean coverage of each 500bp window
    """
    zcat ${regions_bed} | awk -vFS="\\t" -vOFS="\\t" '\$1 != "chrM" && \$1 != "MT" && \$4 > ${threshold}' > ${prefix}.high_coverage.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.high_coverage.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """
}
