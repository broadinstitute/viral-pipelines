version 1.0

import "../tasks/tasks_16S_amplicon.wdl" as qiime_import 

workflow train_classifier {
    meta {
        description: "User imports OTU database that will be trained on your primer sequences."
        author: "Broad Viral Genomics"
        email: "viral_ngs@broadinstitue.org"
        allowNestedInputs: true 
    }
    input {
       File     otu_ref_db
       File     tax_ref_seqs
       String   f_adapter
       String   r_adapter   
    }

    call qiime_import.train_classifier { 
        input: 
            otu_ref_db = otu_ref,
            tax_ref_seqs = taxonomy_ref,
            f_adapter = forward_adapter,
            r_adapter = reverse_adapter
}
} 
