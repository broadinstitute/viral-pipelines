version 1.0

import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_terra.wdl" as terra
import "../tasks/tasks_utils.wdl" as utils

import "demux_deplete.wdl"


workflow demux_deplete_and_table_insert {
    meta {
        description: "Terra only: Demultiplex data from a sequencing run into outputs organized in per-library-lane and per-sample tables."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File          flowcell_tgz
        String?       read_structure

        String?       instrument_model
        String?       sra_title

        File?         sample_rename_map

        File?         collab_ids_tsv
    }

    call demux_deplete.demux_deplete {
        input:
            flowcell_tgz                    = flowcell_tgz,
            instrument_model_user_specified = instrument_model,
            sra_title                       = sra_title,
            read_structure                  = read_structure,
            sample_rename_map               = sample_rename_map
    }

    call terra.check_terra_env

    if(check_terra_env.is_running_on_terra) {
        call terra.create_or_update_sample_tables {
          input:
            flowcell_run_id     = demux_deplete.run_id,
            workspace_name      = check_terra_env.workspace_name,
            workspace_namespace = check_terra_env.workspace_namespace,
            workspace_bucket    = check_terra_env.workspace_bucket_path,

            raw_reads_unaligned_bams     = demux_deplete.raw_reads_unaligned_bams,
            cleaned_reads_unaligned_bams = demux_deplete.cleaned_reads_unaligned_bams,
            meta_by_filename_json        = demux_deplete.meta_by_filename_json,
            meta_by_sample_json          = demux_deplete.meta_by_sample_json
        }
    }

    output {
        Array[File]   raw_reads_unaligned_bams     = demux_deplete.raw_reads_unaligned_bams
        Array[File]   cleaned_reads_unaligned_bams = demux_deplete.cleaned_reads_unaligned_bams
        Array[File]   cleaned_bams_tiny            = demux_deplete.cleaned_bams_tiny
        File          meta_by_filename_json        = demux_deplete.meta_by_filename_json
        Array[Int]    read_counts_raw              = demux_deplete.read_counts_raw
        Array[Int]    read_counts_depleted         = demux_deplete.read_counts_depleted
        File?         sra_metadata                 = demux_deplete.sra_metadata
        File?         cleaned_bam_uris             = demux_deplete.cleaned_bam_uris
        Array[File]   demux_metrics                = demux_deplete.demux_metrics
        Array[File]   demux_commonBarcodes         = demux_deplete.demux_commonBarcodes
        Array[File]   demux_outlierBarcodes        = demux_deplete.demux_outlierBarcodes
        File          multiqc_report_raw           = demux_deplete.multiqc_report_raw
        File          multiqc_report_cleaned       = demux_deplete.multiqc_report_cleaned
        File          spikein_counts               = demux_deplete.spikein_counts
        Int           num_read_files               = length(demux_deplete.cleaned_reads_unaligned_bams)
        String        run_date                     = demux_deplete.run_date
        String        run_id                       = demux_deplete.run_id
    }
}