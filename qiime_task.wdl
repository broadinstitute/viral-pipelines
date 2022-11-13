version 1.0 
##FN 
#QIIME WORKLOW PART 1.A
#1.samtools (universal)
#2. cutadapt (qiime)
#3. vsearch (qiime)

#Part 1A | Step 1: BAM -> QZA
#qiime import bam-(fastq)-> qza 
task qiime_import_from_bam {
    meta{
        description: "Alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format."
    }
    input { 
    File    bam_file
    String  samplename = basename(bam_file, '.bam')
    String  docker = "quay.io/qiime2/core:2022.8" 
    }
    parameter_meta{
        bam_file:   {description: "Input BAM file"}
    }
command <<<
    set -ex -o pipefail 
    #Part 1A | BAM -> FASTQ [Simple samtools command]
    /opt/conda/envs/qiime2-2022.8/bin/samtools fastq -1 $HOME/R1.fastq.gz -2 $HOME/R2.fastq.gz -0 /dev/null ~{bam_file}
    #making new bash variable | regex: (_) -> (-)
    NEWSAMPLENAME=$(echo "~{samplename}" | perl -lape 's/[_ ]/-/g')
    #All names added to one giant file 
    echo ${NEWSAMPLENAME} > NEWSAMPLENAME.txt
    #Make a manifest.txt that contains [1.sample-id 2.R1_fastq 3.R2.fastq]
    #> =overwrite or writes new file 
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > $HOME/manifest.tsv
    #>>= appends 
    #\t= tabs next value 
    echo -e "$NEWSAMPLENAME\t$HOME/R1.fastq.gz\t$HOME/R2.fastq.gz" >> $HOME/manifest.tsv
    #fastq -> bam (provided by qiime tools import fxn)
    /opt/conda/envs/qiime2-2022.8/bin/qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path $HOME/manifest.tsv \
        --input-format PairedEndFastqManifestPhred33V2 \
        --output-path "~{samplename}.qza" 

>>>
output {
    File outfile_qza = "~{samplename}.qza"
    String samplename_master_sheet = read_string("NEWSAMPLENAME.txt")
}

runtime {
    docker: "~{docker}"
    memory: "1GB"
    cpu: 1 
    }
}
#--------------------------------------------
#Part 1 | Step 2:cutadapt: Trim sequences 
#trimreads_trim
#trim = default 
task trimreads{
    meta{
        description:"Removing adapter sequences, primers, and other unwanted sequence from sequence data."
    }
    input {
    File    reads_qza
    String  qza_basename = basename(reads_qza, '.qza')
    #Boolean not_default = false
    String  forward_adapter = "CTGCTGCCTCCCGTAGGAGT"
    String  reverse_adapter = "AGAGTTTGATCCTGGCTCAG"
    Int? min_length    
    String  docker = "quay.io/qiime2/core:2022.8" 
    }
command <<<
    set -ex -o pipefail
    #cutadapt_trim_command
    /opt/conda/envs/qiime2-2022.8/bin/qiime cutadapt trim-paired \
    --i-demultiplexed-sequences ~{reads_qza} \
    --p-front-f ~{forward_adapter} \ 
    --p-front-r ~{reverse_adapter} \ 
    ~{"--p-minimum-length" + min_length} \
    --p-discard-untrimmed \
    --o-trimmed-sequences "~{qza_basename}_trimmed.qza" 
    #trim_visual 
    /opt/conda/envs/qiime2-2022.8/bin/qiime demux summarize \ 
    --i-data "~{qza_basename}_trimmed.qza" \ 
    --o-visualization "~{qza_basename}_trim_summary.qzv"   
>>>
output {
    #trimmed_sequences = paired ends for vsearch
    File    trimmed_sequence_qza = "~{qza_basename}_trimmed.qza"
    File    trimmed_visualization = "~{qza_basename}_trim_summary.qzv" 
}
runtime {
    docker: "~{docker}"
    memory: "1 GB"
    cpu: 1
}
}
#trimreads_untrim
task trimreads_keep_untrim{
    meta{
        description:"Removing adapter sequences, primers, and other unwanted sequence from sequence data."
    }
    input {
    File    reads_qza
    String  qza_basename = basename(reads_qza, '.qza')
    String  forward_adapter = "CTGCTGCCTCCCGTAGGAGT"
    String  reverse_adapter = "AGAGTTTGATCCTGGCTCAG"
    Int? min_length    
    String  docker = "quay.io/qiime2/core:2022.8" 
    }
command <<<
    set -ex -o pipefail
    #cutadapt_trim_command
    /opt/conda/envs/qiime2-2022.8/bin/qiime cutadapt trim-paired \
    --i-demultiplexed-sequences ~{reads_qza} \
    --p-front-f ~{forward_adapter} \ 
    --p-front-r ~{reverse_adapter} \ 
    ~{"--p-minimum-length" + min_length} \
    --p-no-discard-untrimmed \
    --o-trimmed-sequences $HOME/"~{qza_basename}_all_reads.qza" 
    #untrim_visual
    /opt/conda/envs/qiime2-2022.8/bin/qiime demux summarize \ 
    --i-data "~{qza_basename}_all_reads.qza" \ 
    --o-visualization "~{qza_basename}_all_reads_summary.qzv"      

>>>
output {
    #trimmed_sequences = paired ends for vsearch
    File    all_reads_sequences = "~{qza_basename}_all_reads.qza"
    File    all_reads_visualization = "~{qza_basename}_all_reads_summary.qzv"  
}
runtime {
    docker: "~{docker}"
    memory: "1 GB"
    cpu: 1
}
}
#Part 1 | Step 3:VSEARCH: Merge sequences 
task merge_paired_ends { 
    meta{
        description: "Join paired-end sequence reads using vseach's merge_pairs function."
    }
input {
    #Input File: Merge paired reads
    File    trimmed_sequences
    String  trimmed_sequences_basename = basename(trimmed_sequences, '.qza')

    String  docker = "quay.io/qiime2/core:2022.8"
}
command <<< 
    set -ex -o pipefail

    /opt/conda/envs/qiime2-2022.8/bin/qiime vsearch join-paired_end_seqs \
    --i-demultiplexed-seqs ~{trimmed_sequences}.qza \ 
    --o-joined-sequences ~{trimmed_sequences_basename}_joined.qza
    #Visualize 
    /opt/conda/envs/qiime2-2022.8/bin/qiime demux summarize \
    --i-data demux ~{trimmed_sequences}.qza \
    --o-visualization ~{trimmed_sequences_basename}_visualization.qzv
>>>
output {
    File    joined_end_outfile = "~{trimmed_sequences_basename}_joined.qza"
    File    joined_end_visualization = "~{trimmed_sequences_basename}_visualization.qzv"
}
runtime {
    docker: "~{docker}"
    memory: "1 GB"
    cpu: 1
   
}
#Final output: trimmed (or untrimmed depending on user) + merged ends in qza format.
}