version 1.1

import "../tasks/tasks_read_utils.wdl" as tasks_read_utils

workflow fastq_to_ubam {
    meta {
        description: "Convert reads from fastq format (single or paired) to unaligned BAM format."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    call tasks_read_utils.FastqToUBAM
    output {
        File unmapped_bam = FastqToUBAM.unmapped_bam
    }
}
