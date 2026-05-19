version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_centrifuger {
    meta {
        description: "Centrifuger taxonomic classification of a single BAM file. Emits a per-read classification TSV and a Kraken2-style hierarchical kreport."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    call metagenomics.centrifuger

    output {
        File        classification_tsv    = centrifuger.classification_tsv
        File        kreport               = centrifuger.kreport
        Array[File] unclassified_reads    = centrifuger.unclassified_reads
        Array[File] classified_reads      = centrifuger.classified_reads
        String      viralngs_version      = centrifuger.viralngs_version
    }
}
