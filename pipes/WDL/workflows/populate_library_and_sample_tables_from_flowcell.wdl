version 1.0

import "../tasks/tasks_terra.wdl" as terra

workflow populate_library_and_sample_tables_from_flowcell {
    meta {
        description: "Terra only: Populate per-library-lane and per-sample tables from existing demultiplexed flowcell output"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        String flowcell_run_id
    }

    # obtain runtime workspace info necessary to read or change data in
    # Terra tables of the workspace associated with a job
    call terra.check_terra_env

    if(check_terra_env.is_running_on_terra) {
        call terra.create_or_update_sample_tables {
          input:
            flowcell_run_id     = flowcell_run_id,
            workspace_name      = check_terra_env.workspace_name,
            workspace_namespace = check_terra_env.workspace_namespace,
            workspace_bucket    = check_terra_env.workspace_bucket_path
        }
    }
}