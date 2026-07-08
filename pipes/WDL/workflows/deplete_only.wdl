version 1.1

import "../tasks/tasks_taxon_filter.wdl" as taxon_filter

workflow deplete_only {
    meta {
        description: "Taxonomic depletion of reads matching unwanted taxa (such as human)."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call taxon_filter.deplete_taxa

    output {
        File   cleaned_bam               = deplete_taxa.cleaned_bam
        Int    depletion_read_count_pre  = deplete_taxa.depletion_read_count_pre
        Int    depletion_read_count_post = deplete_taxa.depletion_read_count_post
        String viral_classify_version    = deplete_taxa.viralngs_version
    }
}
