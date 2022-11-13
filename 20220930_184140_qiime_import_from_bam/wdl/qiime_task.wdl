version 1.0


task test {
 input { String hello }
command {
     echo "${hello}"
     qiime
}
runtime {
    docker: "quay.io/qiime2/core:2021.2"}
}

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
    String  samplename 
    String  docker = "quay.io/qiime2/core:2021.2" 
    }
    parameter_meta{
        bam_file:   {description: "Input BAM file"}
    }
command <<<
    set -ex -o pipefail 
    #Part 1A | BAM -> FASTQ [Simple samtools command]
    samtools fastq -1 R1.fastq.gz -2 R2.fastq.gz -0 /dev/null ~{bam_file}
    #making new bash variable | regex: (_) -> (-)
    NEWSAMPLENAME=$(echo "~{samplename}" | perl -lape 's/[_ ]/-/g')
    #All names added to one giant file 
    echo $NEWSAMPLENAME > NEWSAMPLENAME.txt
    #Make a manifest.txt that contains [1.sample-id 2.R1_fastq 3.R2.fastq]
    #> =overwrite or writes new file 
    echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest.txt
    #>>= appends 
    #\t= tabs next value 
    echo -e "$NEWSAMPLENAME\tR1.fastq.gz\tR2.fastq.gz" >> manifest.txt
    #fastq -> bam
    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path manifest.txt \
        --output-path "~{samplename}.qza" \
        --input-format PairedEndFastqManifestPhred33V2
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
#Part 1 | Step 2:cutadapt: Trim sequences 
task trimreads{
    meta{
        description:"Removing adapter sequences, primers, and other unwanted sequence from sequence data."
    }
    input {
    #Input File: .qza
    File    reads_qza
    #Boolean not_default = false
    String  forward_adapter = "CTGCTGCCTCCCGTAGGAGT"
    String  reverse_adapter = "AGAGTTTGATCCTGGCTCAG"
    Int? min_length
    Boolean dont_discard_untrimmed = false

    String  docker = "quay.io/qiime2/core:2021.2" 

    }
command <<<
    set -ex -o pipefail
    qiime cutadapt trim-paired \
    --i-demultiplexed-sequences ~{reads_qza}.qza \
    --p-front-f ~{forward_adapter}  \ 
    --p-front-r ~{reverse_adapter} \
    #how to define optional parameter in command block 
    ~{"--p-minimum-length" + min_length} \
    ~{true="--p-discard-untrimmed" false="" dont_discard_untrimmed}
    qiime demux summarize \ 
    --i-data "~{reads_qza}.qza" \ 
    --o-visualization "~{reads_qza}.qzv" 
>>>
output {
    File    seqs_outfile = "~{reads_qza}.qza" 
    File    summary_visulization = "~{reads_qza}.qzv"
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
    File    paired_end_seqs
    String  docker = "quay.io/qiime2/core:2021.2" #check if correct location
}
command <<< 
#VSEARCH 
    #merge paired reads with vsearch
    set -ex -o pipefail

    qiime vsearch join-paired_end_seqs \
    --i-demultiplexed-seqs ~{paired_end_seqs}_paired.qza \ 
    --o-joined-sequences ~{paired_end_seqs}_demux.qza
    #Visualize 
    qiime demux summarize \
    --i-data demux ~{paired_end_seqs}_demux.qza \
    --o-visualization ~{paired_end_seqs}_demux.qzv
>>>
output {
    File    joined_end_outfile = "~{paired_end_seqs}_demux_joined.qza"
    File    joined_end_visualization = "~{paired_end_seqs}_demux_joined.qzv"
}
runtime {
    docker: "~{docker}"
    memory: "1 GB"
    cpu: 1
   
}

}