version 1.0

import "../tasks/tasks_interhost.wdl" as interhost

workflow mafft_and_trim {
    meta {
        description: "MAFFT based multiple alignment followed by trimal-based edge trimming."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call interhost.multi_align_mafft as mafft

    scatter(alignment in mafft.alignments_by_chr) {
        call interhost.trimal_clean_msa as trimal {
            input:
                in_aligned_fasta = alignment
        }
    }

    output {
        Array[File] alignments_by_chr          = mafft.alignments_by_chr
        Array[File] trimmed_alignements_by_chr = trimal.trimal_cleaned_fasta
        File        sampleNamesFile            = mafft.sampleNamesFile
        String      viral_phylo_version        = mafft.viralngs_version
    }
}
