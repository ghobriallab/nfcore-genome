/**
 * Removes reads in high-coverage regions from a BAM/CRAM, mate-aware and in a single pass.
 *
 * The input `keep` BED is the complement of the high-coverage regions (built by
 * BEDTOOLS_KEEP_REGIONS). `samtools view -L keep --fetch-pairs` keeps reads overlapping the
 * kept regions and pulls in their mates, so pairs stay intact (no orphans). Reads in pairs
 * that fall entirely inside a high-coverage region are dropped; a pair straddling a boundary
 * is kept whole. The stream preserves the input coordinate order (no re-sort) and the input
 * format (BAM or CRAM); the output is re-indexed.
 */
process SAMTOOLS_REMOVE_HIGH_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(keep_bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("${prefix}.bam"),  path("${prefix}.bam.bai"),   optional: true, emit: bam
    tuple val(meta), path("${prefix}.cram"), path("${prefix}.cram.crai"), optional: true, emit: cram
    path  "versions.yml",                                                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}.filtered"
    def reference = fasta ? "--reference ${fasta}" : ""
    // Preserve the input format (bam or cram); samtools index picks .crai/.bai accordingly
    file_type     = input.getExtension()
    if ("${input}" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    # Keep reads overlapping the kept (non-high-coverage) regions, fetching mates so pairs
    # stay intact. Single streaming pass; coordinate order is preserved.
    samtools \\
        view \\
        --threads ${task.cpus} \\
        ${reference} \\
        --output-fmt ${file_type} \\
        --fetch-pairs \\
        -L ${keep_bed} \\
        ${args} \\
        -o ${prefix}.${file_type} \\
        ${input}

    samtools index -@ ${task.cpus} ${prefix}.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix    = task.ext.prefix ?: "${meta.id}.filtered"
    file_type = input.getExtension()
    def index_ext = file_type == "cram" ? "crai" : "bai"
    """
    touch ${prefix}.${file_type}
    touch ${prefix}.${file_type}.${index_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
