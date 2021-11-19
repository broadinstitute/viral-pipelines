version 1.0

import "../tasks/tasks_terra.wdl" as terra

workflow update_data_tables {
    meta {
        description: "Create data tables in Terra workspace from provided tsv load file."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call terra.upload_reads_assemblies_entities_tsv

    output {
        Array[String] tables = upload_reads_assemblies_entities_tsv.tables
    }
}