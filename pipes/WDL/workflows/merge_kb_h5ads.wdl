version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow merge_kb_h5ads {
    meta {
        description: "Merges multiple .h5ad files (output from kb_python) into a single .h5ad file."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call metagenomics.kb_merge_h5ads
    output {
        File out_h5ad               = kb_merge_h5ads.kb_merged_h5ad
        String viral_core_version   = kb_merge_h5ads.viralngs_version
    }   
}
