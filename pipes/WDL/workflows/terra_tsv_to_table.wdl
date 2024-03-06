version 1.0

#DX_SKIP_WORKFLOW

import "../tasks/tasks_terra.wdl" as terra

workflow terra_tsv_to_table {
    meta {
        description: "Upload tsv file to Terra data table: insert-or-update on existing rows/columns"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call terra.upload_entities_tsv
}