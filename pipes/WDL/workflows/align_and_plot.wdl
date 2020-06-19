version 1.0

import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_assembly.wdl" as assembly

workflow align_and_plot {
    input {
        String? aligner_options = "-r Random -l 30 -g 40 -x 20 -t 502"
    }

    call assembly.align_reads as align {
        input:
            aligner_options = aligner_options
    }
    call reports.plot_coverage {
        input:
            aligned_reads_bam = align.aligned_only_reads_bam,
            sample_name = basename(basename(basename(align.aligned_only_reads_bam, ".bam"), ".mapped"), ".clean")
    }

    output {
        File   aligned_bam                   = align.aligned_bam
        File   aligned_bam_idx               = align.aligned_bam_idx
        File   aligned_bam_flagstat          = align.aligned_bam_flagstat
        File   aligned_only_reads_bam        = align.aligned_only_reads_bam
        File   aligned_only_reads_bam_idx    = align.aligned_only_reads_bam_idx
        File   aligned_only_reads_fastqc     = align.aligned_only_reads_fastqc
        File   aligned_only_reads_fastqc_zip = align.aligned_only_reads_fastqc_zip
        Int    reads_provided                = align.reads_provided
        Int    reads_aligned                 = align.reads_aligned
        Int    read_pairs_aligned            = align.read_pairs_aligned
        Float  bases_aligned                 = align.bases_aligned
        Float  mean_coverage                 = align.mean_coverage
        String align_viral_core_version      = align.viralngs_version
        File   coverage_plot                 = plot_coverage.coverage_plot
        File   coverage_tsv                  = plot_coverage.coverage_tsv
        Int    reference_length              = plot_coverage.assembly_length
        String plot_viral_core_version       = plot_coverage.viralngs_version
    }
}
