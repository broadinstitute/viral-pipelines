version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_virnucpro_multi {
    meta {
        description: "Runs VirNucPro deep learning viral classification on multiple BAM files in parallel using GPU acceleration."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]+ reads_bams

        Int expected_length = 300
    }

    parameter_meta {
        reads_bams: {
          description: "Set of files containing reads to process. Must be unmapped reads, paired-end or single-end.",
          patterns: ["*.bam"]
        }
        expected_length: {
          description: "Expected sequence length in bp. Must be 300 or 500 to match bundled models (300_model.pth, 500_model.pth). VirNucPro silently accepts other values but produces invalid predictions."
        }
    }

    scatter(reads_bam in reads_bams) {
        call metagenomics.classify_virnucpro {
            input:
                reads_bam = reads_bam,
                expected_length = expected_length
        }
    }

    output {
        Array[File] virnucpro_reports = classify_virnucpro.report_tsv
    }
}
