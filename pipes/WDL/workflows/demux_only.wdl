version 1.0

import "../tasks/tasks_demux.wdl" as tasks_demux
import "../tasks/tasks_reports.wdl" as reports

workflow demux_only {
    meta {
        description: "Picard-based demultiplexing and basecalling from a tarball of a raw BCL directory."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call tasks_demux.illumina_demux

    call reports.MultiQC {
        input:
            input_files = illumina_demux.raw_reads_fastqc_zip
    }

    output {
        Array[File] raw_reads_unaligned_bams = illumina_demux.raw_reads_unaligned_bams
        File        demux_metrics            = illumina_demux.metrics
        File        demux_commonBarcodes     = illumina_demux.commonBarcodes
        File        demux_outlierBarcodes    = illumina_demux.outlierBarcodes
        File        multiqc_report_raw       = MultiQC.multiqc_report
        String      demux_viral_core_version = illumina_demux.viralngs_version
    }
}
