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
        Array[File] input_archive_files
        String?     bucket_path_prefix
        String?     out_dir_name
        
        String?     gcloud_access_token
    }

    parameter_meta {
        input_archive_files: {
            description: "List of input archive files to unpack.",
            patterns: ["*.tar", "*.tar.gz", "*.tgz", "*.tar.bz2", "*.tbz2", "*.tar.xz", "*.txz", "*.tar.lz4", "*.tar.zst"]
        }
        bucket_path_prefix: {
            description: "Path within the Google Storage bucket to unpack the archive files. If not provided, the root of the bucket will be used."
        }
        out_dir_name: {
            description: "Name of the (sub-)directory to unpack the archive contents to within the bucket prefix specified. If not provided, the contents will be unpacked to the bucket prefix."
        }
        gcloud_access_token: {
            description: "Access token for the Google Cloud Storage bucket, for account authorized to write to the bucket specified by 'bucket_path_prefix'. If not provided, the gcloud auth configuration of the execution environment will be obtained via 'gcloud auth print-access-token' for the active authenticated user (on Terra, the service worker/'pet' account)."
        }
    }

    call tasks_terra.check_terra_env

    # only run the task if we are running on GCP or the user provides an auth token to interact with GCP
    # if needed, we can also inspect 'check_terra_env.is_running_on_terra'
    if( check_terra_env.is_backed_by_gcp || defined(gcloud_access_token) ) {
        call tasks_utils.unpack_archive_to_bucket_path {
            input:
                input_archive_files = input_archive_files,
                gcloud_access_token = gcloud_access_token,
                bucket_path_prefix  = if (check_terra_env.is_running_on_terra && check_terra_env.is_backed_by_gcp) then select_first([bucket_path_prefix,check_terra_env.workspace_bucket_path]) else select_first([bucket_path_prefix]),
                out_dir_name        = out_dir_name
        }
    }
}
