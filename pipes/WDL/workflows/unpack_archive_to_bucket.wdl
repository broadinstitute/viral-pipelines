version 1.0

import "../tasks/tasks_utils.wdl" as tasks_utils
import "../tasks/tasks_terra.wdl" as tasks_terra

workflow unpack_archive_to_bucket {
    meta {
        description: "Unpack archive(s) to a target location within a Google Storage bucket"
        author:      "Broad Viral Genomics"
        email:       "viral-ngs@broadinstitute.org"

        allowNestedInputs: true
    }

    input {
        String? gcloud_auth_token
    }

    call tasks_terra.check_terra_env

    if( (check_terra_env.is_running_on_terra && check_terra_env.is_backed_by_gcp) || defined(gcloud_auth_token) ) {
        call tasks_utils.unpack_archive_to_bucket_path
    }
}
