import "tasks_reports.wdl" as reports

workflow align_and_plot_bulk {

    Array[File]+ reads_unmapped_bam_files
    File assembly_fasta
    File? novocraft_license
    
    scatter(reads_unmapped_bam in reads_unmapped_bam_files) {
        call reports.plot_coverage {
            input:
                assembly_fasta = assembly_fasta,
                reads_unmapped_bam = reads_unmapped_bam,
                novocraft_license = novocraft_license,
                
                aligner = "novoalign",
                aligner_options = "-r Random -l 30 -g 40 -x 20 -t 502",
                skip_mark_dupes = false,
                plot_only_non_duplicates = false,
                bin_large_plots = false,
                binning_summary_statistic = "max",
                
                docker = "quay.io/broadinstitute/viral-core"
        }
    }
}
