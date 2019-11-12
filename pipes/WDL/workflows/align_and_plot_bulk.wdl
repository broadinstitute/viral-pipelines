import "tasks_reports.wdl" as reports

workflow align_and_plot_bulk {

    Array[File]+ reads_unmapped_bam_files
    File assembly_fasta
    File? novocraft_license
    
    String?  aligner="novoalign" # novoalign or bwa
    String?  aligner_options="-r Random -l 30 -g 40 -x 20 -t 502"

    Boolean? skip_mark_dupes=false
    Boolean? plot_only_non_duplicates=false
    Boolean? bin_large_plots=false
    String?  binning_summary_statistic="max" # max or min

    String?  docker="quay.io/broadinstitute/viral-core"
  
    
    scatter(reads_unmapped_bam in reads_unmapped_bam_files) {
        call reports.plot_coverage {
            input:
                assembly_fasta = assembly_fasta,
                reads_unmapped_bam = reads_unmapped_bam,
                novocraft_license = novocraft_license,
                
                aligner = aligner,
                aligner_options = aligner_options,
                skip_mark_dupes = skip_mark_dupes,
                plot_only_non_duplicates = plot_only_non_duplicates,
                bin_large_plots = bin_large_plots,
                binning_summary_statistic = binning_summary_statistic,
                
                docker = docker
        }
    }
}
