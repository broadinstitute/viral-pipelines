##Test 
##09.19.22


import "./tasks_qiime_edit_3.wdl" as qiime 

#Part 1:Importing and utilizing Part 1 (Samtools)
workflow qiime_wfl {
    meta {
        description: "Convert BAM files to QZA QIIME workflow."
        author: "fnegrete"
        email:  "viral_ngs@broadinstitute.org"
        allowNestedInputs: true 
    }

    input {
        File    reads_bam
        File    trained_classifier
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
    #_________________________________________
    call qiime.gen_feature_table {
        input: 
            joined_end_outfile = merge_paired_ends.joined_end_outfile
    }
    #_________________________________________
    call qiime.tax_analysis {
        input:
            trained_classifier = trained_classifier,
            rep_seqs_outfile = gen_feature_table.rep_seqs_outfile,
            rep_table_outfile = gen_feature_table.rep_table_outfile
    }
}
