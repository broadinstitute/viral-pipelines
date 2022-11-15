##Test 
##09.19.22

#call qiime.task
#input 
    #input_task_file = #last_task.outfile

#call qiime.task
#input 
    #input_task_file = input_task file (if not file prev exist or first insert)

# ToDo: change relative import based on evential location of tasks file
import "./tasks_16S_amplicon.wdl" as qiime 

#Part 1:Importing and utilizing Part 1 (Samtools)
workflow qiime_wfl {
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
            trimmed_reads_qza = trim_reads.trimmed_reads_qza
    }
    #_________________________________________
    call qiime.gen_feature_table {
        input: 
            joined_end_outfile = merge_paired_ends.joined_end_outfile
    }
    #_________________________________________
    call qiime.train_classifier {
        input:
            otu_ref = otu_ref,
            taxonomy_ref = taxonomy_ref
    }
    #_________________________________________
    call qiime.tax_analysis {
        input:
            trained_classifier = train_classifier.trained_classifier,
            rep-seqs_outfile = gen_feature_table.rep-seqs_outfile,
            rep-table_outfile = gen_feature_table.rep-table_outfile,
    }
}
