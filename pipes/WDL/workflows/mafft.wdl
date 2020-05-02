version 1.0

import "../tasks/tasks_interhost.wdl" as interhost

workflow mafft {
    call interhost.multi_align_mafft
    output {
        File        sampleNamesFile     = multi_align_mafft.sampleNamesFile
        Array[File] alignments_by_chr   = multi_align_mafft.alignments_by_chr
        String      viral_phylo_version = multi_align_mafft.viralngs_version
    }
}
