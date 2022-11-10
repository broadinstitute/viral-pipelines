##Test 
##09.19.22

# ToDo: change relative import based on evential location of tasks file
import "./tasks_qiime.wdl" as qiime 

#Part 1:Importing and utilizing Part 1 (Samtools)
workflow b {
    meta {
        description: "Convert BAM files to QZA QIIME workflow."
        author: "fnegrete"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true 
    }

    input {
        File    reads_bam
        String  sample_name
        Boolean keep_untrimmed_reads 
    }

    call qiime.qiime_import_from_bam {
        input:
            reads_bam  = reads_bam,
            sample_name = sample_name
    }

    #__________________________________________

    call qiime.trim_reads {
        input:
           reads_qza            = qiime_import_from_bam.reads_qza,
           keep_untrimmed_reads = keep_untrimmed_reads
    }

    #__________________________________________
    call qiime.merge_paired_ends {
        input: 
            reads_qza = trim_reads.trimmed_reads_qza
    }
}
