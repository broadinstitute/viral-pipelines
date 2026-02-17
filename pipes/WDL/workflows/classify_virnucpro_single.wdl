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

        Boolean   parallel=false
        Boolean   persist_model=false
        Boolean   use_gpu=false
        Int       batch_size=256
        Int?      esm_batch_size
        Int?      dnabert_batch_size

        Boolean   resume = false
    }

    call metagenomics.classify_virnucpro {
        input:
            reads_input = reads_input,
            expected_length = expected_length,
            parallel = parallel,
            persist_model = persist_model,
            use_gpu = use_gpu,
            batch_size = batch_size,
            esm_batch_size = esm_batch_size,
            dnabert_batch_size = dnabert_batch_size,
            resume = resume
    }

    output {
        File virnuc_pro_scores = classify_virnucpro.virnuc_pro_scores
    }
}
