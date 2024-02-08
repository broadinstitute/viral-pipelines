version 1.0

task trim_rmdup_subsamp {
    meta {
        description: "Trim reads via trimmomatic, remove duplicate reads, and subsample to a desired read count (default of 100,000), bam in, bam out. "
    }
    input { 
        File inBam
        String bam_basename = basename(inBam, '.bam')                    
        File clipDb
        Int n_reads=10000000
        Int machine_mem_gb = 128
        Int cpu = 16
        Int disk_size_gb = 100 
        String docker ="quay.io/broadinstitute/viral-assemble:2.1.33.0"
    }
    parameter_meta {
        inBam: {
            description: "Input BAM file"
            category: "required"
        }
        clipDb: {
            description: "FASTA file that has a list of sequences to trim from the end of reads. These includes various sequencing adapters and primer sequences that may be on the ends of reads, including those for most of the Illumina kits we use."
            category: "required"
        }
        outBam: {
            description: "Cleaned BAM files (default=outbam.bam)"
            category: "other"
        }
        n_reads: {
            description: "Maximum number of reads set to 10000000 by default."
            category: "required"
        }
    }
    command <<<
        set -ex o pipefail
        assembly.py --version | tee VERSION
        #BAM ->FASTQ-> OutBam? https://github.com/broadinstitute/viral-assemble/blob/80bcc1da5c6a0174362ca9fd8bc0b49ee0b4103b/assembly.py#L91
        assembly.py trim_rmdup_subsamp \
        "~{inBam}" \
        "~{clipDb}" \
        "outBam.bam" \
        ~{'--n_reads=' + n_reads}

        #samtools [OutBam -> FASTA]
        #-f 4 (f = include only) (4 = unmapped reads) https://broadinstitute.github.io/picard/explain-flags.html
        samtools fasta "outBam.bam" > "~{bam_basename}.fasta"
    >>>
output {
    File    trimmed_fasta = "~{bam_basename}.fasta"
}
runtime {
    docker:docker
    memory: machine_mem_gb + "GB"
    cpu: cpu
    disks: "local-disk" + disk_size_gb + "LOCAL"
    dx_instance_type: "n2-highmem-16"
}
}

task lca_megablast {
    meta {
        description: "Runs megablast followed by LCA for taxon identification."
    }
    input {
        File    trimmed_fasta
        File    blast_db_tgz
        String  db_name = "copy"
        File    taxonomy_db_tgz
        String  fasta_basename = basename(trimmed_fasta, ".fasta")
        Int     machine_mem_gb = 500 
        Int     cpu = 16
        Int     disk_size_gb = 300
        String  docker = "quay.io/broadinstitute/viral-classify:2.2.3.0"
    }
    parameter_meta {
        trimmed_fasta: {
            description: "Input sequence FASTA file with clean bam reads."
            category: "required"
        }
        blast_db_tgz: {
            description: "Compressed BLAST database."
            category: 'required'
        }
        db_name: {
            description: "BLAST database name (default = nt)."
            category: "other"
        }
        taxonomy_db_tgz: {
            description: "Compressed taxnonomy dataset."
            category: "required"
        }
    }
    command <<<
    # Make directories
    mkdir -p blastdb results taxdump
    read_utils.py extract_tarball \
      ~{blast_db_tgz} blastdb \
      --loglevel=DEBUG
    
    # Unpack taxonomy.dmp
    read_utils.py extract_tarball \
      ~{taxonomy_db_tgz} taxdump \
      --loglevel=DEBUG

    BLASTDB="2blastdb_nt/"
    # Run megablast against nt
    #miniwdl run worked when the Title DB was same as called under db. Remade DB, make sure to note title of DB. 
    blastn -task megablast -query "~{trimmed_fasta}" -db "2blastdb_nt/2nt" -max_target_seqs 50 -num_threads `nproc` -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" -out "~{fasta_basename}.fasta_megablast_nt.tsv"
    
    # Run LCA
    retrieve_top_blast_hits_LCA_for_each_sequence.pl "~{fasta_basename}.fasta_megablast_nt.tsv" taxdump/nodes.dmp 10 > "~{fasta_basename}.fasta_megablast_nt.tsv_LCA.txt"
    # Done
>>>

output {
    File    LCA_output = "~{fasta_basename}.fasta_megablast_nt.tsv_LCA.txt"
}
runtime {
    docker:docker
    memory: machine_mem_gb + "GB"
    cpu: cpu
    disks: "local-disk" + disk_size_gb + "LOCAL"
    dx_instance_type: "n2-highmem-16"
}
}
#blastn -task megablast -query Hep_WGS19_291.downsampled-100000.fasta -db "blastDB" -max_target_seqs 50 -num_threads `nproc` -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" -out "fasta_megablast_nt.tsv"