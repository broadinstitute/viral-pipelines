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
}
