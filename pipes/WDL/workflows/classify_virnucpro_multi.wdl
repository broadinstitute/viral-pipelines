version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_virnucpro_multi {
    meta {
        description: "Runs VirNucPro deep learning viral classification on multiple unmapped reads files in parallel using GPU acceleration."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]+ reads_input_set

        Int       expected_length = 300
    }

    parameter_meta {
        reads_input_set: {
          description: "Set of files containing reads to process. Must be unmapped reads, paired-end or single-end. Accepts BAM or FASTA formats.",
          patterns: ["*.bam", "*.fasta", "*.fa", "*.fasta.gz", "*.fa.gz"]
        }
        expected_length: {
          description: "Expected sequence length in bp. Must be 300 or 500 to match bundled models (300_model.pth, 500_model.pth). VirNucPro silently accepts other values but produces invalid predictions."
        }
    }

    scatter(reads_input in reads_input_set) {
        call metagenomics.classify_virnucpro {
            input:
                reads_input = reads_input,
                expected_length = expected_length
        }
    }

    output {
        Array[File] virnucpro_scores = classify_virnucpro.virnuc_pro_scores
    }
}
