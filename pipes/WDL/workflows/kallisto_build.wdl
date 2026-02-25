version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow kb_build {
    meta {
        description: "Build a kb/kallisto index"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call metagenomics.build_kallisto_db

    output {
        File    kb_reads_extracted      = build_kallisto_db.kb_index
        String  viral_classify_version  = build_kallisto_db.viralngs_version
    }
}
