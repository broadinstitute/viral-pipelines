version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow filter_classified_bam_to_taxa {
    meta {
        description: "Taxonomic filtration of reads utilizing output from a classifier such as kraken1/2/uniq. Can filter out or filter to a specified taxonomic grouping."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call metagenomics.filter_bam_to_taxa
    output {
        File    bam_filtered_to_taxa                        = filter_bam_to_taxa.bam_filtered_to_taxa
        Int     classified_taxonomic_filter_read_count_pre  = filter_bam_to_taxa.classified_taxonomic_filter_read_count_pre
        Int     classified_taxonomic_filter_read_count_post = filter_bam_to_taxa.classified_taxonomic_filter_read_count_post
        String  viral_classify_version                      = filter_bam_to_taxa.viralngs_version
    }
}
