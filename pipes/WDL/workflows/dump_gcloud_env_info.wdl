version 1.0

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
        Boolean is_running_on_terra = check_terra_env.is_running_on_terra
        Boolean is_backed_by_gcp    = check_terra_env.is_backed_by_gcp

        String  workspace_name      = check_terra_env.workspace_name
        String  workspace_namespace = check_terra_env.workspace_namespace

        String  google_project_id   = check_terra_env.google_project_id

        File    env_info            = check_terra_env.env_info
        File    gcloud_config_info  = check_terra_env.gcloud_config_info
    }
}
