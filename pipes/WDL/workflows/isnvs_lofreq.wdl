version 1.0

import "../tasks/tasks_intrahost.wdl" as intrahost

workflow isnvs_lofreq {
    meta {
        description: "variant calls by LoFreq against reference_fasta"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    parameter_meta {
        aligned_bams: {
            description: "reads aligned to the sequence within reference_fasta.",
            patterns: ["*.bam"]
        }
        reference_fasta: {
            description: "Reference sequence to which reads have be aligned",
            patterns: ["*.fasta"]
        }
    }

    input {
        Array[File]+ aligned_bams
        File         reference_fasta
    }

    scatter(aligned_bam in aligned_bams) {
        call intrahost.lofreq as lofreq {
            input:
                aligned_bam = aligned_bam,
                reference_fasta = reference_fasta
        }
    }

    output {
        Array[File] lofreq_vcfs           = lofreq.report_vcf
        String      lofreq_version        = lofreq.lofreq_version[0]
    }
}
