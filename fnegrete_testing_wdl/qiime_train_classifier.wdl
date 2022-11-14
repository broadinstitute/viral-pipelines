#11.14.22
#FN
#Train a classifier. Default is SILVA DB provided already for. 
import "./tasks_qiime_edit_3.wdl" as qiime 

workflow train_classifier {
    meta{
        description: " Training a classifier"
        author: "fnegrete"
        email: "viral_ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File    otu_ref
        File    taxanomy_ref
        String  forward_adapter
        String  reverse_adapter
    }
    call qiime.train_classifier {
        input: 
            otu_ref = otu_ref, 
            taxanomy_ref = taxanomy_ref,
            forward_adapter = forward_adapter, 
            reverse_adapter = reverse_adapter
    }
}