version 1.0

import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_sarscov2.wdl" as sarscov2
import "../tasks/tasks_terra.wdl" as terra
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_utils.wdl" as utils

import "demux_deplete.wdl"


workflow demux_deplete_and_table_insert {
    meta {
        description: "Demultiplex data from a sequencing run into outputs organized in per-library-lane and per-sample tables."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File          flowcell_tgz
        String?       read_structure

        #Array[File]?  biosample_attributes
        String?       instrument_model
        String?        sra_title

        File?         sample_rename_map

        File?         collab_ids_tsv
    }

    # merge biosample attributes tables
    #call utils.tsv_join as biosample_merge {
    #    input:
    #        input_tsvs   = select_all(flatten([biosample_attributes])),
    #        id_col       = 'accession',
    #        out_basename = "biosample_attributes-merged"
    #}
    #call utils.fetch_col_from_tsv as accessioned_samples {
    #  input:
    #    tsv = biosample_merge.out_tsv,
    #    col = 'sample_name'
    #}

    call demux_deplete.demux_deplete {
        input:
            flowcell_tgz                    = flowcell_tgz,
            instrument_model_user_specified = instrument_model,
            sra_title                       = sra_title,
            read_structure                  = read_structure,
            sample_rename_map               = select_first([sample_rename_map])
    }
    String  flowcell_id = demux_deplete.run_id

    ### gather data by biosample
    #call read_utils.group_bams_by_sample {
    #    input:
    #        bam_filepaths = demux_deplete.cleaned_reads_unaligned_bams
    #}

    call terra.check_terra_env

    call terra.create_or_update_sample_tables {
      input:
        flowcell_run_id     = flowcell_id,
        workspace_name      = check_terra_env.workspace_name,
        workspace_namespace = check_terra_env.workspace_namespace,
        workspace_bucket    = check_terra_env.workspace_bucket_path
    }

    #output {
    #    File stdout_log = create_or_update_sample_tables.stdout_log()
    #    File stderr_log = create_or_update_sample_tables.stderr_log()
    #}
}