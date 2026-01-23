version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_virnucpro_single {
    meta {
        description: "Deep learning viral classification of a single BAM file via VirNucPro using GPU acceleration."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File reads_bam
        Int  expected_length = 300
    }

    call metagenomics.classify_virnucpro {
        input:
            reads_bam = reads_bam,
            expected_length = expected_length
    }

    output {
        File report_tsv = classify_virnucpro.report_tsv
    }
}
