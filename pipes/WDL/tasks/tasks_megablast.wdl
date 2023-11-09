version 1.0

task megablast {
    meta {
        description: "Megablast "
    }
    input {
        File    reads_bam
        File    blast_db_tgz
        File    taxonomy_db_tgz
        String  bam_basename = basename(reads_bam, ".bam")
        Int     machine_mem_gb = 128 
        Int     cpu = 16
        Int     disk_size_gb = 300
        String  docker = "quay.io/broadinstitute/viral-classify:2.1.33.0"
    }
    command <<<
    set -ex o pipefail
    #samtools [BAM -> FASTA]
    SAMPLENAME=$(basename reads_bam .bam | perl -lape 's/[_]/-/g')
    echo $SAMPLENAME
    samtools fasta -F 4 reads_bam > "$SAMPLENAME.fasta"
    # Make directories
    cd~
    mkdir -p blastdb queries fasta results blastdb_custom
    # Move input into queries directory
    mv "$SAMPLENAME.fasta" queries
    #tar -xzvf -C for directory 
    tar -xzvf -C "~{blast_db_tgz}" > blast/blastdb/
    # Run megablast against nt
    ncbi/blast \
    blastn -task megablast -query /blast/queries/"$SAMPLENAME.fasta" -db "nt" -max_target_seqs 50 -num_threads 8 \
    -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" \
    -out /blast/results/"$SAMPLENAME.fasta"_megablast_nt.out
    
    # Download nodes.dmp into taxdump directory, remove this since user input and unpack
    cd ~ 
    mkdir taxdump
    cd taxdump
    tar -xvcf -C "~{taxonomy_db_tgz}"
    cd ~
    # Run LCA
    perl retrieve_top_blast_hits_LCA_for_each_sequence.pl results/sample.fasta_megablast_nt.out taxdump/nodes.dmp 10 > results/sample.fasta_megablast_nt.out_LCA.txt
    ​
    # Done
>>>

output {
    File    LCA_output = "~{bam_basename}_LCA.txt"

runtime {
    docker:docker
    memory: machine_mem_gb + "GB"
    cpu: cpu
    disks: "local-disk" + disk_size_gb + "LOCAL"
    instance_type: "n2-highmem-16"
}
}
}