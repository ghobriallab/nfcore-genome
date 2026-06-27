/**
 * Neutralizes pairing flags for orphaned reads left behind by coverage capping.
 *
 * VARIANTBAM (and any coverage-capping/region-removal step) can drop one mate of a
 * proper pair while the surviving read keeps its paired / proper-pair flags. Downstream
 * blat realignment filtering calls `pysam.AlignmentFile.mate()` for every read whose
 * flag claims a proper pair; that lookup raises `ValueError: mate not found` when the
 * mate record is absent, crashing the job on real WGS data.
 *
 * `filter_bam.py` clears the pairing flags of such orphans (turning them into effectively
 * single-end records, preserving sequence/position for coverage) so `mate()` is never
 * called on a missing mate. It requires a coordinate-sorted input, so we sort first, then
 * filter, then re-index.
 */
process FILTER_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0' :
        'community.wave.seqera.io/library/htslib_pysam_samtools:0958fcf8e91c9212' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), path("${prefix}.bam.bai"), emit: bam
    path "versions.yml",                                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}.filtered"
    if ("${bam}" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    # filter_bam.py needs a coordinate-sorted BAM; sort first.
    samtools sort --threads ${task.cpus} -o sorted.bam ${bam}

    filter_bam.py \\
        --in_bam sorted.bam \\
        --out_bam ${prefix}.bam \\
        ${args}

    samtools index -@ ${task.cpus} ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.filtered"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
