import "../tasks/tasks_reports.wdl" as reports

workflow align_and_plot {
    String  reads_unmapped_bam
    call reports.plot_coverage {
        input:
            reads_unmapped_bam = reads_unmapped_bam,
            aligner_options = "-r Random -l 30 -g 40 -x 20 -t 502",
            sample_name = basename(basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt"), ".clean")
    }
}
