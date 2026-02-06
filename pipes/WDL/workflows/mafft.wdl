version 1.1

import "../tasks/tasks_interhost.wdl" as interhost

workflow mafft {
    meta {
        description: "MAFFT multiple-alignment for a set of possibly multi-segment genomes."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call interhost.multi_align_mafft_ref

    output {
        #File        sampleNamesFile     = multi_align_mafft_ref.sampleNamesFile
        Array[File] alignments_by_chr   = multi_align_mafft_ref.alignments_by_chr
    }
}
