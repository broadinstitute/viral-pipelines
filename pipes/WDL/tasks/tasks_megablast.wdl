version 1.0

task trim_rmdup_subsamp {
    meta {
        description: "Trimming unmapped and unaligned BAM sequences"
    }
    input { 
        File inBam
        String bam_basename = basename(bam, '.bam')                    
        File clipDb
        File outBam
        Int n_reads=100000
        String trim_opts = None
        Int machine_mem_gb = 128
        Int disk_size_gb = 100 
        String docker ="quay.io/broadinstitute/viral-assemble:2.1.33.0"
    }
    command <<<
        set -ex o pipefail
        assembly.py --version | tee VERSION
        #BAM ->FASTQ-> OutBam? https://github.com/broadinstitute/viral-assemble/blob/80bcc1da5c6a0174362ca9fd8bc0b49ee0b4103b/assembly.py#L91
        assembly.py trim_rmdup_subsamp_reads \
        ~{'--inBam' + inBamam} \
        ~{'--clipDb'+ clipDb} \
        ~{'--outBam' + outBam} \
        ~{'n_reads' + n_reads} \
        ~{'trimp_opts' + trim_opts}

    #samtools [BAM -> FASTA]
    #-f 4 (f = include only) (4 = unmapped reads) https://broadinstitute.github.io/picard/explain-flags.html
    samtools fasta -f 4 ~{outBam} > "~{bam_basename}.fasta"
    >>>
output {
    File    cleaned_fasta = "~{bam_basename}.fasta"
}
runtime {
    docker:docker
    memory: machine_mem_gb + "GB"
    cpu: cpu
    disks: "local-disk" + disk_size_gb + "LOCAL"
    instance_type: "n2-highmem-16"
}
}


task megablast {
    meta {
        description: "Megablast "
    }
    input {
        File    cleaned_fasta
        File    blast_db_tgz
        File    taxonomy_db_tgz
        String  fasta_basename = basename(cleaned_fasta, ".fasta")
        Int     machine_mem_gb = 128 
        Int     cpu = 16
        Int     disk_size_gb = 300
        String  docker = "quay.io/broadinstitute/viral-classify:2.1.33.0"
    }
    command <<<
    # Make directories
    cd~ 
    mkdir -p blastdb queries fasta results blastdb_custom
    # Move input into queries directory
    mv ~{cleaned_fasta} queries
    #tar -xzvf -C for directory 
    tar -xzvf ~{blast_db_tgz} -C blast/blastdb/
    # Run megablast against nt
    ncbi/blast \
    blastn -task megablast -query /blast/queries/~{cleaned_fasta} -db ~{blast_db_tgz} -max_target_seqs 50 -num_threads 8 \
    -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" \
    -out /blast/results/"~{fasta_basename}.fasta_megablast_nt.out"
    
    # Download nodes.dmp into taxdump directory, remove this since user input and unpack
    cd ~ 
    mkdir taxdump
    cd taxdump
    tar -xzvf ~{taxonomy_db_tgz} 
    cd ~
    # Run LCA
    perl retrieve_top_blast_hits_LCA_for_each_sequence.pl results/sample.fasta_megablast_nt.out taxdump/nodes.dmp 10 > results/sample.fasta_megablast_nt.out_LCA.txt
    ​
    # Done
>>>

output {
    File    LCA_output = "~{cleaned_fasta}_LCA.txt"
}
runtime {
    docker:docker
    memory: machine_mem_gb + "GB"
    cpu: cpu
    disks: "local-disk" + disk_size_gb + "LOCAL"
    instance_type: "n2-highmem-16"
}
}
