version 1.0

import "../tasks/tasks_16S_amplicon.wdl" as qiime_16S

workflow qiime_import_bam {

meta{
    description: "Importing BAM files into QIIME"
    author: "fnegrete"
    email:  "viral_ngs@broadinstitute.org"
    allowNestedInputs: true 
}
input {
    Array[File]    reads_bam
}

call qiime_16S.qiime_import_from_bam {
    input: 
            reads_bam  = reads_bam
    }
} 