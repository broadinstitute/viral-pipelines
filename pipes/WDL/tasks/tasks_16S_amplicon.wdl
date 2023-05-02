version 1.0 

task qiime_import_from_bam {
    meta {
        description: "Parsing demultiplexed fastq BAM files into qiime readable files."
    }
    input { 
        Array[File] reads_bam
        Int    memory_mb = 7000
        Int     cpu = 5
        Int     disk_size_gb = ceil(2*20) + 5
        String  docker     = "quay.io/broadinstitute/qiime2" 
    }
    parameter_meta {
        reads_bam: {
            description: "Unaligned reads in BAM format, one sample per BAM file.",
            category: "required"
        }
        reads_qza: {
            description: "All unaligned reads in a single QZA (QIIME) file.",
            category: "other"
        }
        samplename_master_sheet: {
            description: "File contains all samples names.",
            category: "other"
        }
    }

    command <<<
        set -ex -o pipefail

        #Part 1A | BAM -> FASTQ [Simple samtools command]
        manifest_TSV=manifest.tsv
        echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest.tsv
        for bam in ~{sep=' ' reads_bam}; do
            #making new bash variable | regex: (_) -> (-)
            NEWSAMPLENAME=$(basename $bam .bam  | perl -lape 's/[_]/-/g')
            echo $NEWSAMPLENAME
            samtools fastq -1 $NEWSAMPLENAME.R1.fastq.gz -2 $NEWSAMPLENAME.R2.fastq.gz -0 /dev/null $bam
            #All names added to one giant file 
            echo $NEWSAMPLENAME >> NEWSAMPLENAME.txt
            #>=replaces
            #>>= appends 
            #\t= tabs next value 
            echo -e "$NEWSAMPLENAME\t$(pwd)/$NEWSAMPLENAME.R1.fastq.gz\t$(pwd)/$NEWSAMPLENAME.R2.fastq.gz"
            echo -e "$NEWSAMPLENAME\t$(pwd)/$NEWSAMPLENAME.R1.fastq.gz\t$(pwd)/$NEWSAMPLENAME.R2.fastq.gz" >> manifest.tsv
        done
        # debug
        cat manifest.tsv
        #fastq -> bam (provided by qiime tools import fxn)
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path manifest.tsv \
            --input-format PairedEndFastqManifestPhred33V2 \
            --output-path "batch.qza"
    >>>

    output {
        File   reads_qza               = "batch.qza"
        String samplename_master_sheet = read_string("NEWSAMPLENAME.txt")
    }
    runtime {
        docker: docker
        memory: "~{memory_mb} MiB"
        cpu: cpu
        disk: disk_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}

#Part 1 | Step 2:cutadapt: Trim sequences 
#trimreads_trim
#trim = default 
task trim_reads {

    meta {
        description:"Removing adapter sequences, primers, and other unwanted sequence from sequence data."
    }

    input {
        File    reads_qza
        String  forward_adapter         = "CTGCTGCCTCCCGTAGGAGT"
        String  reverse_adapter         = "AGAGTTTGATCCTGGCTCAG"
        Int     min_length              = 1
        Boolean keep_untrimmed_reads   = false
        Int     memory_mb = 2000
        Int     cpu = 4
        Int     disk_size_gb = ceil(2*size(reads_qza, "GiB")) + 5
        String  docker          = "quay.io/broadinstitute/qiime2" 
    }
    parameter_meta {
        reads_qza: {
            description: "All unaligned reads in a single QZA (QIIME) file.",
            cateogry: "required"
        }
        forward_adapter: {
            description: "Forward amplicon primer sequence.",
            category: "advanced"
        }
        reverse_adapter: {
            description: "Reverse amplicon primer sequence.",
            cateogry: "advanced"
        }
        min_length: {
            description: "Minimum length of the read, cutadapt will discard anything that is shorter than n bp AFTER trimming.Set to default.",
            category: "other"
        }
        keep_untrimmed_reads: {
            description: "Allows you to choose whether or not to discard untrimmed reads.",
            category: "advanced"
        }
        trimmed_reads_qza: {
            description: "Trimmed reads data file.",
            category: "advanced"
        }
        trimmed_visualization: {
            description: "A diagram that compares your demuxed reads before and after cutting (i.e. length of reads, how many reads were retained).",
            category: "advanced"
        }
    }
    command <<<
        set -ex -o pipefail
        qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "~{reads_qza}" \
        --p-front-f "~{forward_adapter}" \
        --p-front-r "~{reverse_adapter}" \
        ~{"--p-minimum-length " + min_length} \
        ~{true='--p-no-discard-untrimmed' false='--p-discard-untrimmed' keep_untrimmed_reads} \
        --o-trimmed-sequences "trimmed.qza"

        #trim_visual 
        qiime demux summarize \
        --i-data "trimmed.qza" \
        --o-visualization "trim_summary.qzv"
    >>>

    output {
        #trimmed_sequences = paired ends for vsearch
        File trimmed_reads_qza     = "trimmed.qza"
        File trimmed_visualization = "trim_summary.qzv" 
    }

    runtime {
        docker: docker
        memory: "${memory_mb} MiB"
        cpu: cpu
        disk: disk_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}

#Part 1 | Step 3:VSEARCH: Merge sequences 
task join_paired_ends { 
    meta {
        description: "Join paired-end sequence reads using vseach's merge_pairs function."
    }
    input {
        #Input File: Merge paired reads
        File    trimmed_reads_qza
        Int     memory_mb = 2000
        Int     cpu = 1
        Int     disk_size_gb = ceil(2*size(trimmed_reads_qza, "GiB")) + 50
        String  docker = "quay.io/broadinstitute/qiime2"
    }
    parameter_meta{
        trimmed_reads_qza: {
            description:"Trimmed reads data file.",
            category: "required"
        }
        joined_end_reads_qza:{
            description: "Merge paired read file.",
            category: "other"
        }
        joined_end_visualization: {
            description: "This summary is especially useful for assessing the length of linked reads and the quality scores at each sequence base position. ",
            category: "other"
        }
    }
    command <<< 
        set -ex -o pipefail
        qiime vsearch join-pairs \
        --i-demultiplexed-seqs ~{trimmed_reads_qza} \
        --o-joined-sequences "joined.qza"

        qiime demux summarize \
        --i-data "joined.qza" \
        --o-visualization "visualization.qzv"
    >>>
    output {
        File joined_end_reads_qza       = "joined.qza"
        File joined_end_visualization = "visualization.qzv"
    }
    runtime {
        docker: docker
        memory: "${memory_mb} MiB"
        cpu: cpu
        disk: disk_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}

task deblur {

    meta {
        description: "Perform sequence quality control for Illumina data using the Deblur workflow with a 16S reference as a positive filter."
        }
    input {
        File    joined_end_reads_qza
        Int     trim_length_var = 300
        Int     memory_mb = 2000
        Int     cpu = 1
        Int     disk_size_gb = ceil(2*size(joined_end_reads_qza, "GiB")) + 5
        String  docker = "quay.io/broadinstitute/qiime2"
    }
    parameter_meta {
        joined_end_reads_qza: {
            description: "Merge paired read file.",
            category: "required"
        }
        trim_length_var: {
            description: "Length that all seqeuences will be trimmed, and discard any sequences that are not at least this long.",
            category: "advanced"
        }
        representative_seqs_qza: {
            description: "Generate a list of the representative sequences. May be useful to the user if they want to blast these sequences or check for correct trimming.",
            category: "other"
        }
        representative_table_qza: {
            description: "Generate a table of the representaitve sequences.",
            category: "other"
        }
        feature_table: {
            description: "A table that represent the number of of features per sample, the number of samples a given feature is found in.",
            category: "other"
        }
        visualize_stats:{
            description: "Generate visualization of deblur stats.",
            category: "other"
        }
    }
        command <<< 
        set -ex -o pipefail

            qiime deblur denoise-16S \
            --i-demultiplexed-seqs ~{joined_end_reads_qza}\
            ~{"--p-trim-length " + trim_length_var} \
            --p-sample-stats \
            --o-representative-sequences "rep_seqs.qza" \
            --o-table "table.qza" \
            --o-stats "stats.qza"
            
            #Generate feature table- give you the number of features per sample 
            qiime feature-table summarize \
            --i-table  "table.qza" \
            --o-visualization   "table.qzv"
            #Generate visualization of deblur stats
            qiime deblur visualize-stats \
            --i-deblur-stats "stats.qza" \
            --o-visualization "stats.qzv"
        >>>
    output {
        File representative_seqs_qza = "rep_seqs.qza"
        File representative_table_qza = "table.qza"
        File feature_table = "table.qzv"
        File visualize_stats = "stats.qzv"

    }
    runtime {
        docker: docker
        memory: "${memory_mb} GB"
        cpu: cpu
        disk: disk_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}
task train_classifier {
    meta {
        descrription: " Upload a classidier trained to classify v1-2 amplicon sequences"
    }
    input {
        File    otu_ref
        File    taxanomy_ref
        String  forward_adapter
        String  reverse_adapter
        Int    min_length  =   100
        Int    max_length  =   500
        String  otu_basename    =   basename(otu_ref, '.qza')
        Int     memory_mb = 2000
        Int     cpu = 1
        Int     disk_size_gb = ceil(2*size(otu_ref, "GiB")) + 5
        String  docker = "quay.io/broadinstitute/qiime2"
    }
    parameter_meta{
        otu_ref: {
            description: "Operational taxonomic units (OTUs) sequences imported as FASTA file.",
            category:"required"
        }
        taxanomy_ref: {
            description: "Reference taxonomy file.",
            category: "required"
        }
        forward_adapter: {
            description: "The forward primer sequence for the amplicon target.",
            category: "advanced"
        }
        reverse_adapter: {
            description: "The reverse primer sequence for the amplicon target.",
            category:"advanced"
        }
        min_length: {
            description: "Minimum length of amplicon sequences.",
            category: "advanced"
        }
        max_length: {
            description: "Maximum length of amplicon sequences.",
            category:"advanced"
        }
        trained_classifier: {
            description: "Trained taxonomic classifier on target amplicon sequences.",
            category: "other"
        }
    }

    command <<<
     set -ex -o pipefail
        CONDA_ENV_NAME=$(conda info --envs -q | awk -F" " '/qiime.*/{ print $1 }')
        conda activate ${CONDA_ENV_NAME}

        qiime tools import \
        --type 'FeatureData[Sequence]' \
        --input-path ~{otu_ref} \
        --output-path "~{otu_basename}_seqs.qza"

        qiime tools import \
        --type 'FeatureData[Taxonomy]'
        --input-format HeaderlessTSVTaxonomyFormat \
        --input-path ~{taxanomy_ref} \
        --output-path "~{otu_basename}_tax.qza"

        qiime feature-classifier extract-reads\
        --i-sequeunces "~{otu_basename}_seqs.qza"\
        --p-f-primer "~{forward_adapter}" \
        --p-r-primer "~{reverse_adapter}" \
        ~{"--p-min-length " + min_length} \
        ~{"--p-max-length " + max_length} \
        --o-reads "~{otu_basename}_v1-2-ref-seqs.qza"

        qiime feature-classifier fit-classifier-naive-bayes \ 
        --i-reference-reads "~{otu_basename}_v1-2-ref-seqs.qza" \ 
        --i-reference-taxonomy "~{otu_basename}_tax.qza" \ 
        --o-classifier "~{otu_basename}_v1-2-classifier.qza"
        >>>
    output {
        File    trained_classifier = "~{otu_basename}_v1-2-classifier.qza"
    }
    runtime {
        docker: docker
        memory: "${memory_mb} GB"
        cpu: cpu
        disk: disk_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}
task tax_analysis {
    meta {
        description: "Protocol describes performing a taxonomic classification with a naive bayes classifier that has been trained on the V1-2 regions amplified by our primers."
    }
    input {
        File    trained_classifier
        File    representative_seqs_qza
        File    representative_table_qza 
        Int     memory_mb = 5
        Int     cpu = 1
        Int     disk_size_gb = 375
        String  docker = "quay.io/broadinstitute/qiime2"
    }
    parameter_meta{ 
        trained_classifier: {
            description: "Trained taxonomic classifier on target amplicon sequences.",
            category: "required"
            }
        representative_seqs_qza: {
            description: "List of representative sequences.",
            category:"required"
            }
        representative_table_qza: {
            description: "Table of representative sequences.",
            category:"other"
            }
        rep_seq_list: {
            description: "Generate list of representative sequences.",
            category:"other"
            }
        tax_classification_graph: {
            description: "Create a bar graph of your taxonomic classification.",
            category:"other"
            }
    }
    command <<<
        set -ex -o pipefail
        qiime feature-classifier classify-sklearn \
        --i-classifier ~{trained_classifier} \
        --i-reads ~{representative_seqs_qza} \
        --o-classification "taxonomy.qza"
        
        qiime feature-table tabulate-seqs \
        --i-data ~{representative_seqs_qza} \
        --o-visualization "list_rep_seqs.qzv"

        qiime taxa barplot \
        --i-table ~{representative_table_qza} \
        --i-taxonomy "taxonomy.qza" \
        --o-visualization "taxa_bar_plots.qzv"
    >>>
    output {
        File rep_seq_list = "list_rep_seqs.qzv" 
        File tax_classification_graph = "taxa_bar_plots.qzv"
}
    runtime {
        docker: docker
        memory: "10 GB"
        cpu: cpu
        disk: disk_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
} 