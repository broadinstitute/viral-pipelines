version 1.0

import "../tasks/tasks_terra.wdl" as terra

workflow terra_table_to_tsv {
    meta {
        description: "Download data table in Terra workspace to tsv file."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call terra.download_entities_tsv

    output {
        File tsv_file = download_entities_tsv.tsv_file
    }
}