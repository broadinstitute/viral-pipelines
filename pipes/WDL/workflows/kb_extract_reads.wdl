version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow kb_extract_reads {
    meta {
        description: "Extracts reads that produced a hit against a kb/kallisto index"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call metagenomics.kb_extract

    output {
        File    kb_extracted_reads       = kb_extract.kb_extract_tar
        String  viral_classify_version   = kb_extract.viralngs_version
    }
}
