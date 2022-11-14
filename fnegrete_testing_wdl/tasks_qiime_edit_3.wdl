version 1.0 
##FN 
#qiime2022.8 version
#QIIME WORKLOW PART 1.A
#1.samtools (universal)
#2. cutadapt 
#3. vsearch 
#4. deblur
#5. 

#Part 1A | Step 1: BAM _> QZA
#qiime import bam_(fastq)_> qza 
task qiime_import_from_bam {

    meta{
        description: "Alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per_position format."
    }

    input { 
        File    reads_bam
        String  sample_name = basename(reads_bam, '.bam')

        String  docker     = "quay.io/qiime2/core:2022.8" 
    }

    parameter_meta {
        reads_bam: {description: "Input BAM file"}
    }

    command <<<
        set _ex _o pipefail
        # obtain the name of the qiime conda environment in the container; flexible about version
        CONDA_ENV_NAME=$(conda info __envs _q | awk _F" " '/qiime.*/{ print $1 }')
        # activate the qiime conda environment
        # seemingly necessary because of:
        #   https://github.com/chanzuckerberg/miniwdl/issues/603
        conda activate ${CONDA_ENV_NAME}

        #Part 1A | BAM _> FASTQ [Simple samtools command]
        samtools fastq _1 $(pwd)/R1.fastq.gz _2 $(pwd)/R2.fastq.gz _0 /dev/null ~{reads_bam}
        #making new bash variable | regex: (_) _> (_)
        NEWSAMPLENAME=$(echo "~{sample_name}" | perl _lape 's/[_]/_/g')
        #All names added to one giant file 
        echo ${NEWSAMPLENAME} > NEWSAMPLENAME.txt
        #Make a manifest.txt that contains [1.sample_id 2.R1_fastq 3.R2.fastq]
        #> =overwrite or writes new file 
        echo _e "sample_id\tforward_absolute_filepath\treverse_absolute_filepath" > manifest.tsv
        #>>= appends 
        #\t= tabs next value 
        echo _e "$NEWSAMPLENAME\t$(pwd)/R1.fastq.gz\t$(pwd)/R2.fastq.gz" >> manifest.tsv
        
        #fastq _> bam (provided by qiime tools import fxn)
        qiime tools import \
            __type 'SampleData[PairedEndSequencesWithQuality]' \
            __input_path manifest.tsv \
            __input_format PairedEndFastqManifestPhred33V2 \
            __output_path "~{sample_name}.qza"
    >>>

    output {
        File   reads_qza               = "~{sample_name}.qza"
        String samplename_master_sheet = read_string("NEWSAMPLENAME.txt")
    }

    runtime {
        docker: "~{docker}"
        memory: "1GB"
        cpu: 1 
    }
}

#____________________________________________
#Part 1 | Step 2:cutadapt: Trim sequences 
#trimreads_trim
#trim = default 
task trim_reads {

    meta {
        description:"Removing adapter sequences, primers, and other unwanted sequence from sequence data."
    }

    input {
        File    reads_qza
        
        String  qza_basename = basename(reads_qza, '.qza')
        #Boolean not_default = false
        String  forward_adapter         = "CTGCTGCCTCCCGTAGGAGT"
        String  reverse_adapter         = "AGAGTTTGATCCTGGCTCAG"
        Int    min_length              = 1
        Boolean? keep_untrimmed_reads   = false

        String  docker          = "quay.io/qiime2/core:2022.8" 
    }

    command <<<
        set _ex _o pipefail
        # obtain the name of the qiime conda environment in the container; flexible about version
        CONDA_ENV_NAME=$(conda info __envs _q | awk _F" " '/qiime.*/{ print $1 }')
        # activate the qiime conda environment
        # seemingly necessary because of:
        #   https://github.com/chanzuckerberg/miniwdl/issues/603
        conda activate ${CONDA_ENV_NAME}

        qiime cutadapt trim_paired \
        __i_demultiplexed_sequences "~{reads_qza}" \
        __p_front_f "~{forward_adapter}" \
        __p_front_r "~{reverse_adapter}" \
        ~{"__p_minimum_length " + min_length} \
        ~{true='__p_no_discard_untrimmed' false='__p_discard_untrimmed' keep_untrimmed_reads} \
        __o_trimmed_sequences "~{qza_basename}_trimmed.qza"

        #trim_visual 
        qiime demux summarize \
        __i_data "~{qza_basename}_trimmed.qza" \
        __o_visualization "~{qza_basename}_trim_summary.qzv"
    >>>

    output {
        #trimmed_sequences = paired ends for vsearch
        File trimmed_reads_qza     = "~{qza_basename}_trimmed.qza"
        File trimmed_visualization = "~{qza_basename}_trim_summary.qzv" 
    }

    runtime {
        docker: "~{docker}"
        memory: "10 GB"
        cpu: 1
    }
}

#Part 1 | Step 3:VSEARCH: Merge sequences 
task merge_paired_ends { 
    meta {
        description: "Join paired_end sequence reads using vseach's merge_pairs function. Perform sequence quality control for Illumina data using the Deblur workflow with a 16S reference as a positive filter."
        }
    input {
        #Input File: Merge paired reads
        File    trimmed_reads_qza
        String  reads_basename = basename(trimmed_reads_qza, '.qza')

        String  docker = "quay.io/qiime2/core:2022.8"
    }

    command <<< 
        set _ex _o pipefail
        # obtain the name of the qiime conda environment in the container; flexible about version
        CONDA_ENV_NAME=$(conda info __envs _q | awk _F" " '/qiime.*/{ print $1 }')
        # activate the qiime conda environment
        # seemingly necessary because of:
        #   https://github.com/chanzuckerberg/miniwdl/issues/603
        conda activate ${CONDA_ENV_NAME}

        qiime vsearch join_pairs \
        __i_demultiplexed_seqs "~{reads_basename}_trimmed.qza" \
        __o_joined_sequences "~{reads_basename}_joined.qza"

        qiime demux summarize \
        __i_data "~{reads_basename}_joined.qza" \
        __o_visualization "~{reads_basename}_visualization.qzv"
    >>>
    output {
        File joined_end_outfile       = "~{reads_basename}_joined.qza"
        File joined_end_visualization = "~{reads_basename}_visualization.qzv"
    }
    runtime {
        docker: "~{docker}"
        memory: "10 GB"
        cpu: 1
    }
#Final output: trimmed (or untrimmed depending on user) + merged ends in qza format.
}

task gen_feature_table {

    meta {
        description: "Perform sequence quality control for Illumina data using the Deblur workflow with a 16S reference as a positive filter."
        }
    input {
        File    joined_end_outfile
        String  joined_end_basename = basename(joined_end_outfile, '.qza')
        Int    trim_length_var = 300     
        String  docker = "quay.io/qiime2/core:2022.8"
    }
    command <<< 
            set _ex _o pipefail
            CONDA_ENV_NAME=$(conda info __envs _q | awk _F" " '/qiime.*/{ print $1 }')
            conda activate ${CONDA_ENV_NAME}
            ##activated conda environments 
            qiime deblur denoise_16S \
            __i_demultiplexed_seqs "~{joined_end_basename}_demux_joined.qza" \ 
            ~{"__p_trim_length" + trim_length_var} \ 
            __p_sample_stats \
            __o_representative_sequences "~{joined_end_basename}_rep_seqs.qza" \
            __o_table "~{joined_end_basename}_table.qza"
            
            #Generate feature table_ give you the number of features per sample 
            qiime feature_table summarize \
            __i_table  "~{joined_end_basename}_table.qza" \
            __o_visualization   "~{joined_end_basename}_table.qzv"
            #Generate visualization of deblur stats
            qiime deblur visualize_stats \
            __i_deblur_stats "~{joined_end_basename}_deblur_stats.qza" \
            __o_visualization "~{joined_end_basename}_deblur_stats.qzv"
        >>>
    output {
        #how many output files do i need
        File rep_seqs_outfile = "~{joined_end_basename}_rep_seqs.qza"
        File rep_table_outfile = "~{joined_end_basename}_table.qza"
        File feature_table = "~{joined_end_basename}_table.qzv"
        File visualize_stats = "~{joined_end_basename}_deblur_stats.qzv"

    }
    runtime {
        docker: "~{docker}"
        memory: "10 GB"
        cpu: 1
    } 
}
task train_classifier {
    meta {
        descrription: " Upload a classidier trained to classify v1_2 amplicon sequences"
    }
    input {
        File    otu_ref
        File    taxanomy_ref
        String  forward_adapter
        String  reverse_adapter
        Int    min_length  =   100
        Int    max_length  =   500
        String  otu_basename    =   basename(otu_ref, '.qza')
        String  docker = "quminiay.io/qiime2/core:2022.8"
    }
    command <<<
     set _ex _o pipefail
        CONDA_ENV_NAME=$(conda info __envs _q | awk _F" " '/qiime.*/{ print $1 }')
        conda activate ${CONDA_ENV_NAME}
        ##activated conda environments 
        
        #is this otu different than the one above? or is it the rep_seqs outifle?
        #is this how you specify a path?
        qiime tools import \
        __type 'FeatureData[Sequence]' \
        __input_path ~{otu_ref} \
        __output_path "~{otu_basename}_seqs.qza"

        #is this db different than the one above?
        qiime tools import \
        __type 'FeatureData[Taxonomy]'
        __input_format HeaderlessTSVTaxonomyFormat \
        __input_path ~{taxanomy_ref} \
        __output_path "~{otu_basename}_tax.qza"

        qiime feature_classifier extract_reads\
        __i_sequeunces "~{otu_basename}_seqs.qza"\
        __p_f_primer "~{forward_adapter}"\
        __p_r_primer "~{reverse_adapter}"\
        ~{"__p_min_length" + min_length}\
        ~{"__p_max_length" + max_length}\
        __o_reads "~{otu_basename}_v1_2_ref_seqs.qza"
#might have to be broken down into two tasks
        qiime feature_classifier fit_classifier_naive_bayes\ 
        __i_reference_reads "~{otu_basename}_v1_2_ref_seqs.qza"\ 
        __i_reference_taxonomy "~{otu_basename}_tax.qza"\ 
        __o_classifier "~{otu_basename}_v1_2_classifier.qza"
    >>>
    output {
        File    trained_classifier = "~{otu_basename}_v1_2_classifier.qza"
    }
    runtime {
        docker: "~{docker}"
        memory: "10 GB"
        cpu: 1
    }

}
task tax_analysis {
    meta {
        description: "Protocol describes performing a taxonomic classification with a naive bayes classifier that has been trained on the V1_2 regions amplified by our primers."
    }
    input {
        File trained_classifier
        File rep_seqs_outfile
        File rep_table_outfile 
        String  basename  =   basename(trained_classifier, '.qza')
        String  docker = "quay.io/qiime2/core:2022.8"
    }
    command <<<
    set _ex _o pipefail
        CONDA_ENV_NAME=$(conda info __envs _q | awk _F" " '/qiime.*/{ print $1 }')
        conda activate ${CONDA_ENV_NAME}
        ##activated conda environments 
        
        qiime feature_classifier classify_sklearn \
        __i_classifier ~{trained_classifier}\
        __i_reads ~{rep_seqs_outfile}\
        __o_classifcation "~{basename}_tax.qza"
        
        qiime feature_table tabulate_seqs \
        __i_date ~{rep_seqs_outfile}\
        __o_visualization "~{basename}_rep_seqs.qzv"

        qiime taxa barplot\
        __i_table ~{rep_table_outfile}\
        __i_taxonomy "~{basename}_tax.qza"\
        __o_visualization "~{basename}_bar_plots.qzv"
    >>>
    output {
        File rep_seq_list = "~{basename}_rep_seqs.qzv" 
        File tax_classification_graph = "~{basename}_bar_plots.qzv"
    }
    runtime {
        docker: "~{docker}"
        memory: "10 GB"
        cpu: 1
    }
}
# One caviat is that you need to visualize the bar plot using qiime2.org 