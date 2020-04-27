version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow filter_classified_bam_to_taxa {
    call metagenomics.filter_bam_to_taxa
}
