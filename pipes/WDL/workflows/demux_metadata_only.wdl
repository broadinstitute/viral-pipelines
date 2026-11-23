version 1.1

#DX_SKIP_WORKFLOW

import "../tasks/tasks_demux.wdl" as demux
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_terra.wdl" as terra
import "../tasks/tasks_utils.wdl" as utils

workflow demux_metadata_only {
    meta {
        description: "Picard-based demultiplexing and basecalling from a tarball of a raw BCL directory, followed by QC metrics, depletion, and SRA submission prep."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]+        samplesheets
        File?               sample_rename_map

        Array[File]         raw_reads_unaligned_bams
        Array[File]         cleaned_reads_unaligned_bams

        Map[String,String]  run_info
        File                meta_by_filename_json
        Array[File]         biosample_map_tsvs
        String?             instrument_model_user_specified
        String              sra_title

        Boolean             insert_demux_outputs_into_terra_tables=false
    }

    scatter(samplesheet in samplesheets) {
        call demux.samplesheet_rename_ids {
            input:
                old_sheet  = samplesheet,
                rename_map = sample_rename_map
        }
    }

    if(insert_demux_outputs_into_terra_tables) {
        call terra.check_terra_env

        if(check_terra_env.is_running_on_terra) {
            call terra.create_or_update_sample_tables {
              input:
                flowcell_run_id     = run_info['run_id'],
                workspace_name      = check_terra_env.workspace_name,
                workspace_namespace = check_terra_env.workspace_namespace,

                raw_reads_unaligned_bams     = raw_reads_unaligned_bams,
                cleaned_reads_unaligned_bams = cleaned_reads_unaligned_bams,
                meta_by_filename_json        = meta_by_filename_json
            }
        }
    }

    if(length(biosample_map_tsvs) > 0) {
        # NCBI biosample metadata is available

        #### merge biosample attribute tsvs (iff provided with more than one)
        if (length(biosample_map_tsvs) > 1) {
            call utils.tsv_join as biosample_map_tsv_join {
                input:
                    input_tsvs   = biosample_map_tsvs,
                    id_col       = 'accession',
                    out_suffix   = ".tsv",
                    out_basename = "biosample-attributes-merged"
            }
        }
        File biosample_map_tsv = select_first(flatten([[biosample_map_tsv_join.out_tsv], biosample_map_tsvs]))

        #### biosample metadata mapping
        call ncbi.biosample_to_table {
            input:
                biosample_attributes_tsv = biosample_map_tsv,
                raw_bam_filepaths        = raw_reads_unaligned_bams,
                demux_meta_json          = meta_by_filename_json
        }

        #### SRA submission prep
        call ncbi.sra_meta_prep {
            input:
                cleaned_bam_filepaths = cleaned_reads_unaligned_bams,
                biosample_map         = biosample_map_tsv,
                library_metadata      = samplesheet_rename_ids.new_sheet,
                platform              = "ILLUMINA",
                paired                = (run_info['indexes'] == '2'),

                out_name              = "sra_metadata-~{run_info['run_id']}.tsv",
                instrument_model      = select_first(flatten([[instrument_model_user_specified],[run_info['sequencer_model']]])),
                title                 = sra_title
        }

        if(insert_demux_outputs_into_terra_tables && select_first([check_terra_env.is_running_on_terra])) {
            call terra.upload_entities_tsv as terra_load_biosample_data {
                input:
                    workspace_name   = select_first([check_terra_env.workspace_name]),
                    terra_project    = select_first([check_terra_env.workspace_namespace]),
                    tsv_file         = biosample_to_table.sample_meta_tsv
            }
        }
    }


    output {
        File?       sra_metadata                             = sra_meta_prep.sra_metadata
        
        String      instrument_model_inferred                = select_first(flatten([[instrument_model_user_specified],[run_info['sequencer_model']]]))

        File?       terra_library_table                      = create_or_update_sample_tables.library_metadata_tsv
        File?       terra_sample_library_map                 = create_or_update_sample_tables.sample_membership_tsv
        File?       terra_sample_metadata                    = biosample_to_table.sample_meta_tsv
    }
}
