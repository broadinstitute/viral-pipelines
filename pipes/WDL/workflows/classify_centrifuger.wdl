version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_centrifuger {
    meta {
        description: "Centrifuger taxonomic classification of a single BAM file. Emits a per-read classification TSV, a Kraken2-style hierarchical kreport, and a Kraken2-compatible per-read report only when k == 1."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    call metagenomics.centrifuger

    if (length(centrifuger.kraken2_reads_reports) > 0) {
        File extracted_kraken2_reads_report = centrifuger.kraken2_reads_reports[0]
    }

    output {
        File        classification_tsv    = centrifuger.classification_tsv
        File        kreport               = centrifuger.kreport
        File?       kraken2_reads_report  = extracted_kraken2_reads_report
        String      viralngs_version      = centrifuger.viralngs_version
    }
}
