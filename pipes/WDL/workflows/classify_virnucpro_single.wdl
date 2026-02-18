version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_virnucpro_single {
    meta {
        description: "Deep learning viral classification of a single unmapped reads file via VirNucPro using GPU acceleration."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File reads_input
        Int  expected_length = 300
    }

    call metagenomics.classify_virnucpro {
        input:
            reads_input = reads_input,
            expected_length = expected_length
    }

    output {
        File virnuc_pro_scores = classify_virnucpro.virnuc_pro_scores
    }
}
