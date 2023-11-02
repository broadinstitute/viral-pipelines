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
        String  docker = "insert_docker_image"
    }
    command <<<
    set -ex o pipefail
    #samtools [BAM -> FASTA]
    SAMPLENAME=$(basename reads_bam .bam | perl -lape 's/[_]/-/g')
    echo $SAMPLENAME
    samtools fasta -F 4 reads_bam > "$SAMPLENAME.fasta"
    # Make directories
    cd ; mkdir -p blastdb queries fasta results blastdb_custom
    # Move input into queries directory
    mv "$SAMPLENAME.fasta" queries
    # [REMOVED] Download nt (15 minutes) - should be already added from the Docker?
    # Run megablast against nt
    docker run \
    -v $HOME/blastdb:/blast/blastdb:ro -v $HOME/blastdb_custom:/blast/blastdb_custom:ro \
    -v $HOME/queries:/blast/queries:ro \
    -v $HOME/results:/blast/results:rw \
    ncbi/blast \
    blastn -task megablast -query /blast/queries/sample.fasta -db "nt" -max_target_seqs 50 -num_threads 8 \
    -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" \
    -out /blast/results/sample.fasta_megablast_nt.out
    
    # Download nodes.dmp into taxdump directory
    cd ; mkdir taxdump ; cd taxdump
    curl ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz > taxdump.tar.gz
    tar -xf taxdump.tar.gz
    cd
    ​
    # Run LCA
    perl retrieve_top_blast_hits_LCA_for_each_sequence.pl results/sample.fasta_megablast_nt.out taxdump/nodes.dmp 10 > results/sample.fasta_megablast_nt.out_LCA.txt
    ​
    # Done
>>>

output {
    File    LCA_output = "~{bam_basename}_LCA.txt"
    File    sample_fasta = "sample.fasta"
}
runtime {
    docker:docker
    memory: machine_mem_gb + "GB"
    cpu: cpu
    disks: "local-disk" + disk_size_gb + "LOCAL"
    instance_type: "n2-highmem-16"
#do we want to add max retries here?
}
}