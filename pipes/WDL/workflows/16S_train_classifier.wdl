version 1.0

import "../tasks/tasks_16S_amplicon.wdl" as qiime 

workflow train_classifier {
    meta {
        description: "User imports OTU database that will be trained on your primer sequences."
        author: "Broad Viral Genomics"
        email: "viral_ngs@broadinstitue.org"
        allowNestedInputs: true 
    }
    input {
       File     otu_ref
       File     taxonomy_ref
       String   forward_adapter
       String   reverse_adapter   
    }

    call qiime.train_classifier { 
        input: 
            otu_ref = otu_ref,
            taxonomy_ref = taxonomy_ref,
            forward_adapter = forward_adapter,
            reverse_adapter = reverse_adapter
}
} 
