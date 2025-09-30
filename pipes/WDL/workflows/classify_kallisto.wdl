version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_kallisto {
    meta {
        description: "Taxonomic classification of RNA-seq data via kallisto."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call metagenomics.kallisto

    output {
        File    kallisto_count_tar = kallisto.kallisto_count_tar
        String  viral_classify_version = kallisto.viralngs_version
    }
}
