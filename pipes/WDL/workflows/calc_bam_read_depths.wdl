version 1.0

import "../tasks/tasks_read_utils.wdl" as read_utils

workflow calc_bam_read_depths {
    meta {
        description: "Generates read depth tables."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    parameter_meta {
        aligned_bams: {
            description: "Reads aligned to a reference sequence.",
            patterns: ["*.bam"]
        }
    }

    input {
        Array[File]+ aligned_bams
    }

    scatter(aligned_bam in aligned_bams) {
        
        call read_utils.read_depths as depth {
            input:
                aligned_bam = aligned_bam
        }
    }

    output {
        Array[File] read_depths = depth.read_depths
    }
}
