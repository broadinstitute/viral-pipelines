version 1.1

#DX_SKIP_WORKFLOW

import "../tasks/tasks_terra.wdl" as terra

workflow dump_gcloud_env_info {
    meta {
        description: "Write system and gcloud environment info to output files."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call terra.check_terra_env

    output {
        Boolean is_running_on_terra     = check_terra_env.is_running_on_terra
        Boolean is_backed_by_gcp        = check_terra_env.is_backed_by_gcp

        String  google_project_id       = check_terra_env.google_project_id

        String  workspace_uuid          = check_terra_env.workspace_uuid
        String  workspace_name          = check_terra_env.workspace_name
        String  workspace_namespace     = check_terra_env.workspace_namespace
        String  workspace_bucket_path   = check_terra_env.workspace_bucket_path

        String  method_version          = check_terra_env.method_version
        String  method_source           = check_terra_env.method_source
        String  method_path             = check_terra_env.method_path
        
        String  input_table_name        = check_terra_env.input_table_name
        String  input_row_id            = check_terra_env.input_row_id

        String  top_level_submission_id = check_terra_env.top_level_submission_id

        File    env_info                = check_terra_env.env_info
        File    gcloud_config_info      = check_terra_env.gcloud_config_info
    }
}
