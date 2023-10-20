version 1.0

#DX_SKIP_WORKFLOW

import "../tasks/tasks_terra.wdl" as terra

workflow dump_gcloud_env_info {
    meta {
        description: "Write system and gcloud environment info to output files."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call terra.get_gcloud_env_info

    output {
        Array[String] gcloud_env_info_files     = get_gcloud_env_info.env_info_files
        File          additional_command_stdout = get_gcloud_env_info.additional_command_stdout
    }
}
