version 1.0

import "../tasks/tasks_16S_amplicon.wdl" as qiime

workflow train_classifier_16S {
    meta {
        description: "User imports OTU database that will be trained on your primer sequences."
        author: "Broad Viral Genomics"
        email: "viral-ngs@broadinstitue.org"
        allowNestedInputs: true 
    }
    input {
       File     otu_ref
       File     taxanomy_ref
       String   forward_adapter
       String   reverse_adapter   
    }

    call qiime.train_classifier { 
        input: 
            otu_ref = otu_ref,
            taxanomy_ref = taxanomy_ref,
            forward_adapter = forward_adapter,
            reverse_adapter = reverse_adapter
    }
} 
