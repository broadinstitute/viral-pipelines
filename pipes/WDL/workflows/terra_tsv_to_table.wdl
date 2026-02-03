version 1.1

#DX_SKIP_WORKFLOW

import "../tasks/tasks_terra.wdl" as terra
import "../tasks/tasks_utils.wdl" as utils

workflow terra_tsv_to_table {
    meta {
        description: "Upload tsv files to Terra data table: insert-or-update on existing rows/columns"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File?]+ tsv_files
    }

    call terra.check_terra_env

    call utils.cat_except_headers {
        input:
            infiles = select_all(tsv_files),
            out_filename = "terra_upload.tsv"
    }

    call terra.upload_entities_tsv {
        input:
            tsv_file         = cat_except_headers.out_tsv,
            workspace_name   = check_terra_env.workspace_name,
            terra_project    = check_terra_env.workspace_namespace
    }
}
