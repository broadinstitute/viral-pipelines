version 1.0 

task qiime_import_from_bam {
    meta {
        description: "Parsing demultiplexed fastq BAM files into qiime readable files."
    }
    input { 
        File    reads_bam
        String  sample_name
        Int     memory_mb = 2000
        Int     cpu = 1
        Int     disk_size_gb = ceil(2*size(reads_bam, "GiB")) + 5
        String  docker     = "quay.io/qiime2/core:2022.8" 
    }
    parameter_meta {
        reads_bam: {description: "Input BAM file"}
    }

    command <<<
        set -ex -o pipefail
        #testing to see if source is correct
        # obtain the name of the qiime conda environment in the container flexible about version
        #https://github.com/chanzuckerberg/miniwdl/issues/603

        # DEBUGGING: print the name of the current shell (bash, sh, zsh, csh, etc.)
        CURRENT_SHELL_NAME="$(ps -hp $$ | grep -v 'PID' | awk '{print $4}' | sed 's/\-//g')"
        echo ${CURRENT_SHELL_NAME}

        # DEBUGGING: print whether the current shell is an interactive (i.e. login) session 
        [ -z "$PS1" ] && echo "This shell is not interactive" || echo "This shell is interactive"

        # DEBUGGING: print environment variables
        env

        # DEBUGGING: print the name of the current user
        whoami

        # add the (base) miniconda bin directory to the system PATH list and re-index executables
        export PATH="/opt/conda/conda/bin:${PATH}"; hash -r

        # call conda init to ensure conda startup calls are in the user's shell rc dotfile
        # 
        conda init ${CURRENT_SHELL_NAME} -v

        # start a new [login] shell (to re-load shell rc dotfile)
        exec "$CURRENT_SHELL_NAME" -l

        # DEBUGGING: === now in a new shell ===

        # DEBUGGING: print the name of the current shell (bash, sh, zsh, csh, etc.)
        ps -hp $$ | grep -v 'PID' | awk '{print $4}' | sed ('s/\-//g')
        
        # DEBUGGING: print whether the current shell is an interactive (i.e. login) session 
        [ -z "$PS1" ] && echo "This shell is not interactive" || echo "This shell is interactive"

        # DEBUGGING: print environment variables
        env

        # add the (base) miniconda bin directory to the system PATH list and re-index executables
        # needed a second time since this is a new shell after the call 'exec "$CURRENT_SHELL_NAME" -l'
        export PATH="/opt/conda/conda/bin:${PATH}"; hash -r

        # determine name of the qiime conda environment available in this container (variable with version)
        QIIME_CONDA_ENV_NAME=$(conda info --envs -q | awk -F" " '/qiime.*/{ print $1 }')

        # print the qiime conda environment name
        echo "${QIIME_CONDA_ENV_NAME}"

        # activate the qiime conda environment
        conda activate ${QIIME_CONDA_ENV_NAME}

        # DEBUGGING: print environment variables
        env

        #Part 1A | BAM -> FASTQ [Simple samtools command]
        samtools fastq -1 $(pwd)/R1.fastq.gz -2 $(pwd)/R2.fastq.gz -0 /dev/null ~{reads_bam}
        #making new bash variable | regex: (_) -> (-)
        NEWSAMPLENAME=$(echo "~{sample_name}" | perl -lape 's/[_]/-/g')
        #All names added to one giant file 
        echo ${NEWSAMPLENAME} > NEWSAMPLENAME.txt
        #Make a manifest.txt that contains [1.sampleid 2.R1_fastq 3.R2.fastq]
        #> =overwrite or writes new file 
        echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest.tsv
        #>>= appends 
        #\t= tabs next value 
        echo -e "$NEWSAMPLENAME\t$(pwd)/R1.fastq.gz\t$(pwd)/R2.fastq.gz" >> manifest.tsv
        
        #fastq -> bam (provided by qiime tools import fxn)
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path manifest.tsv \
            --input-format PairedEndFastqManifestPhred33V2 \
            --output-path "~{sample_name}.qza"
    >>>

    output {
        File   reads_qza               = "~{sample_name}.qza"
        String samplename_master_sheet = read_string("NEWSAMPLENAME.txt")
    }
    runtime {
        docker: docker
        memory: "${memory_mb} MiB"
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
        
        String  qza_basename = basename(reads_qza, '.qza')
        #Boolean not_default = false
        String  forward_adapter         = "CTGCTGCCTCCCGTAGGAGT"
        String  reverse_adapter         = "AGAGTTTGATCCTGGCTCAG"
        Int    min_length              = 1
        Boolean keep_untrimmed_reads   = false
        Int     memory_mb = 2000
        Int     cpu = 1
        Int     disk_size_gb = ceil(2*size(reads_qza, "GiB")) + 5
        String  docker          = "quay.io/qiime2/core:2022.8" 
    }

    command <<<
        set -ex -o pipefail
        # obtain the name of the qiime conda environment in the container; flexible about version
        CONDA_ENV_NAME=$(conda info --envs -q | awk -F" " '/qiime.*/{ print $1 }')
        # activate the qiime conda environment
        # seemingly necessary because of:
        # https://github.com/chanzuckerberg/miniwdl/issues/603
        conda activate ${CONDA_ENV_NAME}

        qiime cutadapt trim-paired \
        --i-demultiplexed-sequences "~{reads_qza}" \
        --p-front-f "~{forward_adapter}" \
        --p-front-r "~{reverse_adapter}" \
        ~{"--p-minimum-length " + min_length} \
        ~{true='--p-no-discard-untrimmed' false='--p-discard-untrimmed' keep_untrimmed_reads} \
        --o-trimmed-sequences "~{qza_basename}_trimmed.qza"

        #trim_visual 
        qiime demux summarize \
        --i-data "~{qza_basename}_trimmed.qza" \
        --o-visualization "~{qza_basename}_trim_summary.qzv"
    >>>

    output {
        #trimmed_sequences = paired ends for vsearch
        File trimmed_reads_qza     = "~{qza_basename}_trimmed.qza"
        File trimmed_visualization = "~{qza_basename}_trim_summary.qzv" 
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
        String  reads_basename = basename(trimmed_reads_qza, '.qza')
        Int     memory_mb = 2000
        Int     cpu = 1
        Int     disk_size_gb = ceil(2*size(trimmed_reads_qza, "GiB")) + 5
        String  docker = "quay.io/qiime2/core:2022.8"
    }

    command <<< 
        set -ex -o pipefail
        # obtain the name of the qiime conda environment in the container; flexible about version
        CONDA_ENV_NAME=$(conda info --envs -q | awk -F" " '/qiime.*/{ print $1 }')
        # activate the qiime conda environment
        # seemingly necessary because of:
        # https://github.com/chanzuckerberg/miniwdl/issues/603
        conda activate ${CONDA_ENV_NAME}

        qiime vsearch join-pairs \
        --i-demultiplexed-seqs ~{trimmed_reads_qza} \
        --o-joined-sequences "~{reads_basename}_joined.qza"

        qiime demux summarize \
        --i-data "~{reads_basename}_joined.qza" \
        --o-visualization "~{reads_basename}_visualization.qzv"
    >>>
    output {
        File joined_end_reads_qza       = "~{reads_basename}_joined.qza"
        File joined_end_visualization = "~{reads_basename}_visualization.qzv"
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
        String  joined_end_basename = basename(joined_end_reads_qza, '.qza')
        Int     trim_length_var = 300
        Int     memory_mb = 2000
        Int     cpu = 1
        Int     disk_size_gb = ceil(2*size(joined_end_reads_qza, "GiB")) + 5
        String  docker = "quay.io/qiime2/core:2022.8"
    }
        command <<< 
        set -ex -o pipefail
        # obtain the name of the qiime conda environment in the container; flexible about version
        CONDA_ENV_NAME=$(conda info --envs -q | awk -F" " '/qiime.*/{ print $1 }')
        # activate the qiime conda environment
        # seemingly necessary because of:
        # https://github.com/chanzuckerberg/miniwdl/issues/603
        conda activate ${CONDA_ENV_NAME}
            qiime deblur denoise-16S \
            --i-demultiplexed-seqs ~{joined_end_reads_qza} \
            ~{"--p-trim-length " + trim_length_var} \
            --p-sample-stats \
            --o-representative-sequences "~{joined_end_basename}_rep_seqs.qza" \
            --o-table "~{joined_end_basename}_table.qza" \
            --o-stats "~{joined_end_basename}_stats.qza"
            
            #Generate feature table- give you the number of features per sample 
            qiime feature-table summarize \
            --i-table  "~{joined_end_basename}_table.qza" \
            --o-visualization   "~{joined_end_basename}_table.qzv"
            #Generate visualization of deblur stats
            qiime deblur visualize-stats \
            --i-deblur-stats "~{joined_end_basename}_stats.qza" \
            --o-visualization "~{joined_end_basename}_stats.qzv"
        >>>
    output {
        File representative_seqs_qza = "~{joined_end_basename}_rep_seqs.qza"
        File representative_table_qza = "~{joined_end_basename}_table.qza"
        File feature_table = "~{joined_end_basename}_table.qzv"
        File visualize_stats = "~{joined_end_basename}_stats.qzv"

    }
    runtime {
        docker: docker
        memory: "${memory_mb} MiB"
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
        String  docker = "quay.io/qiime2/core:2022.8"
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
        memory: "${memory_mb} MiB"
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
        String  basename  =   basename(trained_classifier, '.qza')
        Int     memory_mb = 2000
        Int     cpu = 1
        Int     disk_size_gb = ceil(2*size(trained_classifier, "GiB")) + 5
        String  docker = "quay.io/qiime2/core:2022.8"
    }
    command <<<
        set -ex -o pipefail
        # obtain the name of the qiime conda environment in the container; flexible about version
        CONDA_ENV_NAME=$(conda info --envs -q | awk -F" " '/qiime.*/{ print $1 }')
        # activate the qiime conda environment
        # seemingly necessary because of:
        #   https://github.com/chanzuckerberg/miniwdl/issues/603
        conda activate ${CONDA_ENV_NAME}
        
        qiime feature-classifier classify-sklearn \
        --i-classifier ~{trained_classifier} \
        --i-reads ~{representative_seqs_qza} \
        --o-classification "~{basename}_tax.qza"
        
        qiime feature-table tabulate-seqs \
        --i-data ~{representative_seqs_qza} \
        --o-visualization "~{basename}_rep_seqs.qzv"

        qiime taxa barplot \
        --i-table ~{representative_table_qza} \
        --i-taxonomy "~{basename}_tax.qza" \
        --o-visualization "~{basename}_bar_plots.qzv"
    >>>
    output {
        File rep_seq_list = "~{basename}_rep_seqs.qzv" 
        File tax_classification_graph = "~{basename}_bar_plots.qzv"
}
    runtime {
        docker: docker
        memory: "${memory_mb} MiB"
        cpu: cpu
        disk: disk_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}