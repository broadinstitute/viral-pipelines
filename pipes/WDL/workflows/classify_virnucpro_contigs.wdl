version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_virnucpro_contigs {
    meta {
        description: "Classify contigs as viral or non-viral based on confidence-weighted delta scores from VirNucPro chunk-level predictions."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File    virnucpro_tgz

        Float   min_viral_prop    = 0.1
        Float   min_nonviral_prop = 0.1
        Int     min_chunks        = 5
        String  id_col            = "Modified_ID"
        String  id_pattern        = "(NODE_\\d+)"
    }

    call metagenomics.classify_virnucpro_contigs as classify_contigs {
        input:
            virnucpro_tgz        = virnucpro_tgz,
            min_viral_prop       = min_viral_prop,
            min_nonviral_prop    = min_nonviral_prop,
            min_chunks           = min_chunks,
            id_col               = id_col,
            id_pattern           = id_pattern
    }

    output {
        File virnucpro_contigs_classified = classify_contigs.contig_classifications
    }
}
