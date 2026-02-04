version 1.0

import "../tasks/tasks_demux.wdl" as tasks_demux
import "../tasks/tasks_reports.wdl" as reports

workflow demux_only {
    meta {
        description: "Picard-based demultiplexing and basecalling from a tarball of a raw BCL directory."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        String? instrument_model_user_specified
    }

    call tasks_demux.illumina_demux

    call reports.multiqc_from_bams as multiqc {
        input:
            input_bams = illumina_demux.raw_reads_unaligned_bams
    }

    output {
        Array[File] raw_reads_unaligned_bams  = illumina_demux.raw_reads_unaligned_bams
        File        demux_metrics             = illumina_demux.metrics
        File        demux_commonBarcodes      = illumina_demux.commonBarcodes
        File        demux_outlierBarcodes     = illumina_demux.outlierBarcodes
        File        multiqc_report_raw        = multiqc.multiqc_report
        Array[File] fastqcs                   = multiqc.fastqc_html

        String             run_date           = illumina_demux.run_info['run_start_date']
        Map[String,String] run_info           = illumina_demux.run_info
        File               run_info_json      = illumina_demux.run_info_json
        String             run_id             = illumina_demux.run_info['run_id']
        String             run_lane_count     = illumina_demux.flowcell_lane_count

        String      instrument_model_inferred = select_first(flatten([[instrument_model_user_specified],[illumina_demux.run_info['sequencer_model']]]))
        String      demux_viral_core_version  = illumina_demux.viralngs_version

        Map[String,Map[String,String]] meta_by_filename      = illumina_demux.meta_by_filename
        Map[String,Map[String,String]] meta_by_sample        = illumina_demux.meta_by_sample
        File                           meta_by_filename_json = illumina_demux.meta_by_filename_json
        File                           meta_by_sample_json   = illumina_demux.meta_by_sample_json
    }
}
