/**
 * Removes reads overlapping high-coverage regions from a CRAM/BAM.
 * `samtools view -L <bed>` selects reads overlapping the regions; `-U` writes the
 * *unselected* reads (everything NOT overlapping) to the output, effectively
 * dropping reads in the high-coverage windows. The kept output is re-indexed and
 * keeps the same format (BAM or CRAM) as the input.
 */
process SAMTOOLS_REMOVE_HIGH_COVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(bed)
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
    # Keep only reads that do NOT overlap the high-coverage regions.
    # -L selects reads overlapping the BED; -U sends the rest (the ones we keep) to the output.
    # The selected reads are discarded to /dev/null (no --write-index here, it would try to
    # index /dev/null). The kept -U output is indexed in a separate step below.
    samtools \\
        view \\
        --threads ${task.cpus} \\
        ${reference} \\
        --output-fmt ${file_type} \\
        -L ${bed} \\
        -U ${prefix}.${file_type} \\
        ${args} \\
        -o /dev/null \\
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
