//
// QC on CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_STATS                } from '../../../modules/nf-core/samtools/stats/main'
include { MOSDEPTH                      } from '../../../modules/nf-core/mosdepth/main'
include { MOSDEPTH_HIGH_COVERAGE        } from '../../../modules/local/mosdepth_high_coverage/main'
include { BEDTOOLS_KEEP_REGIONS         } from '../../../modules/local/bedtools/keep_regions/main'
include { SAMTOOLS_REMOVE_HIGH_COVERAGE } from '../../../modules/local/samtools/remove_high_coverage/main'

workflow CRAM_QC_MOSDEPTH_SAMTOOLS {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai ]
    fasta     // channel: [mandatory] [ fasta ]
    fasta_fai // channel: [mandatory] [ meta, fasta_fai ]
    intervals

    main:
    versions = Channel.empty()
    reports = Channel.empty()

    // Reports run on cram
    SAMTOOLS_STATS(cram, fasta)

    MOSDEPTH(cram.combine(intervals.map { meta, bed -> [bed ?: []] }), fasta)

    // Parse the 500bp-window coverage BED to extract high-coverage regions (mean cov > threshold).
    // Gated off on the post-recalibration call via ext.when in conf/modules/modules.config,
    // so it only runs after the first (post-dedup) mosdepth.
    MOSDEPTH_HIGH_COVERAGE(MOSDEPTH.out.regions_bed)

    // Only act on BEDs that actually contain high-coverage regions (non-empty file).
    // Empty BEDs never launch the removal steps.
    high_cov_bed = MOSDEPTH_HIGH_COVERAGE.out.bed.filter { _meta, bed -> bed.size() > 0 }

    // Build the complement (regions to keep) from the high-coverage BED + genome sizes,
    // then remove reads in the high-coverage regions (mate-aware, single pass).
    BEDTOOLS_KEEP_REGIONS(high_cov_bed, fasta_fai)

    cram_to_filter = cram
        .join(BEDTOOLS_KEEP_REGIONS.out.bed, failOnDuplicate: true)

    SAMTOOLS_REMOVE_HIGH_COVERAGE(cram_to_filter, fasta, fasta_fai)

    // Build the final alignment: swap in the filtered file for samples that had high-coverage
    // regions, keep the original alignment for everyone else. Removal preserves BAM/CRAM format,
    // so mix both output paths (only one is populated per sample).
    filtered = SAMTOOLS_REMOVE_HIGH_COVERAGE.out.bam.mix(SAMTOOLS_REMOVE_HIGH_COVERAGE.out.cram)
    alignment = cram
        .join(filtered, remainder: true)
        .map { row ->
            // matched: [meta, orig_file, orig_index, filt_file, filt_index]; unmatched: [meta, orig_file, orig_index]
            def (meta, orig_file, orig_index, filt_file, filt_index) = row
            filt_file ? [ meta, filt_file, filt_index ] : [ meta, orig_file, orig_index ]
        }

    // Gather all reports generated
    reports = reports.mix(SAMTOOLS_STATS.out.stats)
    reports = reports.mix(MOSDEPTH.out.global_txt)
    reports = reports.mix(MOSDEPTH.out.regions_txt)

    // Gather versions of all tools used
    versions = versions.mix(MOSDEPTH.out.versions)
    versions = versions.mix(SAMTOOLS_STATS.out.versions)
    versions = versions.mix(MOSDEPTH_HIGH_COVERAGE.out.versions)
    versions = versions.mix(BEDTOOLS_KEEP_REGIONS.out.versions)
    versions = versions.mix(SAMTOOLS_REMOVE_HIGH_COVERAGE.out.versions)

    emit:
    reports
    high_coverage_bed = MOSDEPTH_HIGH_COVERAGE.out.bed // channel: [ meta, bed ]
    alignment                                          // channel: [ meta, file, index ] — filtered where high-cov regions existed, original otherwise
    versions                                           // channel: [ versions.yml ]
}
