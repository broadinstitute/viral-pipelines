version 1.0

import "../tasks/tasks_qiime_import_bam.wdl" as import

workflow qiime_import_bam {
    
meta{
    description: "Importing BAM files into QIIME"
    author: "fnegrete"
    email:  "viral_ngs@broadinstitute.org"
    allowNestedInputs: true 
}
input {
    File    reads_bam
    File    sample_name
}

call import.qiime_import_bam {
    input: 
            reads_bam  = reads_bam,
            sample_name = sample_name
}
} 