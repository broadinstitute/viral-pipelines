version 1.0

import "../tasks/tasks_16S_amplicon.wdl" as qiime

workflow train_16S_classifier {
    meta {
        description: "User imports OTU database that will be trained on your primer sequences. All outputs can be visualized using https://view.qiime2.org/"
        author: "Broad Viral Genomics"
        email: "viral-ngs@broadinstitue.org"
        allowNestedInputs: true 
    }
    call qiime.train_classifier{
        input:
            otu_ref = otu_ref,
            taxonomy_ref = taxonomy_ref,
            forward_adapter = forward_adapter, 
            reverse_adapter = reverse_adapter

    }
} 
