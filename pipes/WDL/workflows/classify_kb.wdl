version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_kallisto {
    meta {
        description: "Taxonomic classification of RNA-seq data via kb_python"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call metagenomics.kb

    output {
        File    kb_extract_tar = kb_extract.kb_extract_tar
        String  viral_classify_version = kb_extract.viralngs_version
    }
}
