version 1.0

import "../tasks/tasks_interhost.wdl" as interhost

workflow trimal {
    meta {
        description: "Trim a multiple sequence alignment with Trimal."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call interhost.trimal_clean_msa

    output {
        File trimmed_alignment = trimal_clean_msa.trimal_cleaned_fasta
    }
}
