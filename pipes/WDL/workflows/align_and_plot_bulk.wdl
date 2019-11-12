import "tasks_reports.wdl" as reports

workflow align_and_plot_bulk {

    File assembly_fasta
    File reads_unmapped_bam
    File? novocraft_license
    
    call reports.plot_coverage {
        input:
            aligner_options = "-r Random -l 30 -g 40 -x 20 -t 502",
            assembly_fasta = assembly_fasta,
            reads_unmapped_bam = reads_unmapped_bam,
            novocraft_license = novocraft_license
    }
}
