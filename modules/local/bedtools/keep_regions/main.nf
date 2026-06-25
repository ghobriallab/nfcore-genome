/**
 * Builds the "keep" regions for high-coverage read removal: the complement of the
 * high-coverage BED across the whole genome. The input BED is sorted and merged first
 * (mosdepth --by windows can be unsorted/overlapping/adjacent), then complemented against
 * the genome sizes derived from the .fai so every contig — including those with no
 * high-coverage regions — is represented.
 */
process BEDTOOLS_KEEP_REGIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"

    input:
    tuple val(meta),  path(bed)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("${prefix}.keep.bed"), emit: bed
    path  "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # genome sizes from the .fai (chrom, length), sorted to match the BED sort order
    cut -f1,2 ${fai} | sort -k1,1 > genome.txt

    # sort + merge the high-coverage regions, then complement to get the regions to keep
    sort -k1,1 -k2,2n ${bed} \\
        | bedtools merge -i - \\
        | bedtools complement -i - -g genome.txt \\
        > ${prefix}.keep.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.keep.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
