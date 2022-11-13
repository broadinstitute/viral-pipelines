##Test 
##09.19.22

import "../qiime_task.wdl" as toolbox 
#Part 1:Importing and utilizing Part 1 (Samtools)
workflow b {
    meta {
        description: "Convert BAM files to QZA QIIME workflow."
        author: "fnegrete"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true 
    }
    input {
        File bam_file
        String samplename
        Boolean keep_untrim 
    }
    call toolbox.qiime_import_from_bam {
        input:
            bam_file = bam_file,
            samplename = samplename
    }
    #__________________________________________

    if(!keep_untrim) {
        call toolbox.trimreads
            input:
               reads_qza = qiime_import_from_bam.outfile_qza
        }
    if(keep_untrim) {
        call toolbox.trimreads_keep_untrim
            input:
                reads_qza = qiime_import_from_bam.outfile_qza
        }

    #__________________________________________
   call toolbox.merge_paired_ends {
        input: 
            trimmed_sequences = trimreads.trimmed_sequence_qza
    }
}
