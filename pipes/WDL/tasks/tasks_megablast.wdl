version 1.0

task trim_rmdup_subsamp {
    meta {
        description: "Trim reads via trimmomatic, remove duplicate reads, and subsample to a desired read count (default of 100,000), bam in, bam out. "
    }
    input { 
        File inBam
        String bam_basename = basename(inBam, '.bam')
        File outBam = "outbam.bam"                    
        File clipDb
        Int n_reads=10000000
        Int machine_mem_gb = 128
        Int cpu = 16
        Int disk_size_gb = 100 
        String docker ="quay.io/broadinstitute/viral-assemble:2.3.1.3"
    }
    parameter_meta {
        inBam: {
            description: "Input BAM file",
            category: "required"
        }
        clipDb: {
            description: "FASTA file that has a list of sequences to trim from the end of reads. These includes various sequencing adapters and primer sequences that may be on the ends of reads, including those for most of the Illumina kits we use.",
            category: "required"
        }
        outBam: {
            description: "Cleaned BAM files (default=outbam.bam)",
            category: "other"
        }
        n_reads: {
            description: "Maximum number of reads set to 10000000 by default.",
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
        "~{outBam}" \
        ~{'--n_reads=' + n_reads}

        #samtools [OutBam -> FASTA]
        #-f 4 (f = include only) (4 = unmapped reads) https://broadinstitute.github.io/picard/explain-flags.html
        samtools fasta "~{outBam}" > "~{bam_basename}.fasta"
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
        File    taxdb
        String  db_name
        File    taxonomy_db_tgz
        String  fasta_basename = basename(trimmed_fasta, ".fasta")
        Int     machine_mem_gb = 500 
        Int     cpu = 16
        Int     disk_size_gb = 300
        String  docker = "quay.io/broadinstitute/viral-classify:2.2.4.2"
    }
    parameter_meta {
        trimmed_fasta: {
            description: "Input sequence FASTA file with clean bam reads.",
            category: "required"
        }
        blast_db_tgz: {
            description: "Compressed BLAST database.",
            category: 'required'
        }
        db_name: {
            description: "BLAST database name (default = nt).",
            category: "other"
        }
        taxonomy_db_tgz: {
            description: "Compressed taxnonomy dataset.",
            category: "required"
        }
    }
    command <<<
    #Extract BLAST DB tarball
    read_utils.py extract_tarball \
      ~{blast_db_tgz} . \
      --loglevel=DEBUG
    
    # Extract taxonomy DB tarball 
    read_utils.py extract_tarball \
      ~{taxonomy_db_tgz} . \
      --loglevel=DEBUG

    '''
    #Extract taxid map file tarball
    read_utils.py extract_tarball \
        ~{taxdb} . \
        --loglevel=DEBUG
    '''

    #Set permissions 
    chmod +x /opt/viral-ngs/source/retrieve_top_blast_hits_LCA_for_each_sequence.pl
    chmod +x /opt/viral-ngs/source/LCA_table_to_kraken_output_format.pl
    # Verify BLAST database
    blastdbcmd -db "~{db_name}" -info
    if [ $? -ne 0 ]; then
        echo "Database '~{db_name}' not found or is inaccessible."
        exit 1
    else
        echo "Database '~{db_name}' found and accessible."
    fi
    #miniwdl run worked when the Title DB was same as called under db. Remade DB, make sure to note title of DB. 
    #Log start time 
    START_TIME=$(date +%s)
    # Run megablast against nt
    blastn -task megablast -query "~{trimmed_fasta}" -db "~{db_name}" -max_target_seqs 50 -num_threads `nproc` -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" -out "~{fasta_basename}.fasta_megablast_nt.tsv"
    #Log end time
    END_TIME=$(date +%s)
     # Calculate elapsed time
    ELAPSED_TIME=$(($END_TIME - $START_TIME))
    echo "BLAST step took $ELAPSED_TIME seconds." > blast_elapsed_time.txt
    # Run LCA
    retrieve_top_blast_hits_LCA_for_each_sequence.pl "~{fasta_basename}.fasta_megablast_nt.tsv" nodes.dmp 10 > "~{fasta_basename}.fasta_megablast_nt.tsv_LCA.txt"
    # Run Krona output conversion 
    LCA_table_to_kraken_output_format.pl "~{fasta_basename}.fasta_megablast_nt.tsv_LCA.txt" "~{trimmed_fasta}" > "~{fasta_basename}.kraken.tsv"
    # Done
>>>

output {
    File    LCA_output = "~{fasta_basename}.fasta_megablast_nt.tsv_LCA.txt"
    File    kraken_output_fromat =  "~{fasta_basename}.kraken.tsv"
    File    elapsed_time_normal_blastn = "blast_elapsed_time.txt"
}
runtime {
    docker:docker
    memory: machine_mem_gb + "GB"
    cpu: cpu
    disks: "local-disk" + disk_size_gb + "HDD"
    dx_instance_type: "n2-highmem-16"
}
}

task blastoff {
    meta{
        description:"Blastoff wrapper"
    }
    input{
        File    trimmed_fasta
        Int     host_species = 9606
        Int     stage2_min_id = 98
        Int     stage2_min_qcov = 98 
        #Are these id/qcov b/w stage 2 & 3 ever different?
        Int     stage3_min_id = 98
        Int     stage3_min_qcov = 98
        File    blast_db_tgz
        File    taxonomy_db_tgz
        String  db_name
        String  fasta_basename = basename(trimmed_fasta, ".fasta")
        Int     machine_mem_gb = 500 
        Int     cpu = 2
        Int     disk_size_gb = 300
        String  docker = "quay.io/broadinstitute/viral-classify:2.2.4.2"

    }
    command <<<
    #Extract BLAST DB tarball
    read_utils.py extract_tarball \
      ~{blast_db_tgz} . \
      --loglevel=DEBUG

    # Extract taxonomy DB tarball (includes nodes.dmp)
    read_utils.py extract_tarball \
      ~{taxonomy_db_tgz} . \
      --loglevel=DEBUG
    #permissions 
    chmod +x /opt/viral-ngs/source/retrieve_top_blast_hits_LCA_for_each_sequence.pl
    chmod +x /opt/viral-ngs/source/retrieve_most_common_taxonids_in_LCA_output.pl
    chmod +x /opt/viral-ngs/source/filter_LCA_matches.pl
    chmod +x /opt/viral-ngs/source/add_one_value_column.pl
    chmod +x /opt/viral-ngs/source/retrieve_sequences_appearing_or_not_appearing_in_table.pl
    chmod +x /opt/viral-ngs/source/generate_LCA_table_for_sequences_with_no_matches.pl
    chmod +x /opt/viral-ngs/source/concatenate_tables.pl
    chmod +x /opt/viral-ngs/source/LCA_table_to_kraken_output_format.pl


    #STAGE 1 
    #Subsamples 100 random reads from original FASTA file
    select_random_sequences.pl "~{trimmed_fasta}" 100 > "~{fasta_basename}_subsampled.fasta"
    #run megablast on random reads x nt
    #switched output from out to txt for readability issues
    blastn -task megablast -query "~{fasta_basename}_subsampled.fasta"  -db "~{db_name}" -max_target_seqs 50 -num_threads `nproc` -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" -out "~{fasta_basename}_subsampled.fasta_megablast_nt.tsv" 
    # Run LCA
    retrieve_top_blast_hits_LCA_for_each_sequence.pl "~{fasta_basename}_subsampled.fasta_megablast_nt.tsv" nodes.dmp 1 1 > "~{fasta_basename}_subsampled.fasta_megablast_nt.tsv_LCA.txt"
    #Looks for most frequently matched taxon IDs and outputs a list
    retrieve_most_common_taxonids_in_LCA_output.pl "~{fasta_basename}_subsampled.fasta_megablast_nt.tsv_LCA.txt" species 10 1 > "sample_specific_db_taxa.txt" 
    # Create an empty sample_specific_db_taxa.txt if it doesn't exist
    touch sample_specific_db_taxa.txt
    #adding host_species to sample_specific_db_taxa.txt
    echo "~{host_species}" >> "sample_specific_db_taxa.txt"
    #ensure file is sorted and unique 
    sort sample_specific_db_taxa.txt | uniq > sample_specific_db_taxa_unique.txt
    mv sample_specific_db_taxa_unique.txt sample_specific_db_taxa.txt
    echo "input sequences to stage 2:"
    grep ">" "~{trimmed_fasta}" | wc -l
    echo "--END STAGE 1"

    #STAGE 2
    echo "---START STAGE 2"
    echo "megablast sample-specific database start"
    #Run blastn w/ taxidlist specific 
    blastn -task megablast -query "~{fasta_basename}_subsampled.fasta" -db "~{db_name}" -max_target_seqs 50 -num_threads `nproc` -taxidlist "sample_specific_db_taxa.txt" -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" -out "~{fasta_basename}_megablast_sample_specific_db.tsv" 
    #Run LCA on last output
    retrieve_top_blast_hits_LCA_for_each_sequence.pl  "~{fasta_basename}_megablast_sample_specific_db.tsv" nodes.dmp 2 >  "~{fasta_basename}_megablast_sample_specific_db_LCA.txt" 
    #filter
    filter_LCA_matches.pl "~{fasta_basename}_megablast_sample_specific_db_LCA.txt" 1 0 0 "~{stage2_min_id}" 999 "~{stage2_min_qcov}" 999 > "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt"
    #add one clmn value: database
    add_one_value_column.pl "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt" "database" "sample-specific" > "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt_column_added.txt"
    mv "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt_column_added.txt" "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt"
    #add one clmn value: classified
    add_one_value_column.pl "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt" "classified" "classified" > "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt_column_added.txt"
    mv "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt_column_added.txt" "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt"
    #retrieves collection of unclassified sequences
    retrieve_sequences_appearing_or_not_appearing_in_table.pl "~{trimmed_fasta}" "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt" 0 0 > "~{fasta_basename}_megablast_sample_specific_db_unclassified.fasta"
    # megablast_sample_specific_db_${sample_fasta}_unclassified.fasta = "~{fasta_basename}_megablast_sample_specific_db_unclassified.fasta"
    echo "input sequences to stage 3"
    grep ">" "~{fasta_basename}_megablast_sample_specific_db_unclassified.fasta" | wc -l 
    echo "--END STAGE 2"
    echo "---START STAGE 3"
    
    #Stage 3 
    #/blast/results/${sample_fasta}_megablast_sample_specific_db_megablast_nt.out = "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.tsv"
    blastn -task megablast -query "~{fasta_basename}_megablast_sample_specific_db_unclassified.fasta" -db "~{db_name}" -max_target_seqs 50 -num_threads `nproc` -outfmt "6 qseqid sacc stitle staxids sscinames sskingdoms qlen slen length pident qcovs evalue" -out "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.tsv"
    retrieve_top_blast_hits_LCA_for_each_sequence.pl "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.tsv" nodes.dmp 10 > "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt"
    filter_LCA_matches.pl "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt" 0 0 0 "~{stage3_min_id}" 999 "~{stage3_min_qcov}" 999 > "~{fasta_basename}_megablast_sample_specific_db_megablast_nt_LCA_classified.txt"
    #add one column: database, nt
    add_one_value_column.pl "~{fasta_basename}_megablast_sample_specific_db_megablast_nt_LCA_classified.txt" "database" "nt" > "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt_column_added.txt"
    mv "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt_column_added.txt" "~{fasta_basename}_megablast_sample_specific_db_megablast_nt_LCA_classified.txt"
    #add one column: classified, classified
    add_one_value_column.pl "~{fasta_basename}_megablast_sample_specific_db_megablast_nt_LCA_classified.txt" "classified" "classified" > "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt_column_added.txt"
    mv "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_classified.txt_column_added.txt" "~{fasta_basename}_megablast_sample_specific_db_megablast_nt_LCA_classified.txt"
    filter_LCA_matches.pl "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt" 0 0 0 "~{stage3_min_id}" 999 "~{stage3_min_qcov}" 999 1 > "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt"
    add_one_value_column.pl "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt" "database" "nt" > "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt_column_added.txt"
    mv "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt_column_added.txt" "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt"
    add_one_value_column.pl "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt" "classified" "unclassified" > "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt_column_added.txt"
    mv "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt_column_added.txt" "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt"
    generate_LCA_table_for_sequences_with_no_matches.pl "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt" "~{fasta_basename}_megablast_sample_specific_db_unclassified.fasta" > "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt"
    add_one_value_column.pl "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt" "database" "nt" > "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt_column_added.txt"
    mv "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt_column_added.txt" "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt" 
    add_one_value_column.pl "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt" "classified" "unclassified (no matches)" > "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt_column_added.txt"
    mv "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt_column_added.txt" "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt" 
    #Table 1 (classified sequences from stage 2) "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt"
    #Table 2 (classified sequences from stage 3) "~{fasta_basename}_megablast_sample_specific_db_megablast_nt_LCA_classified.txt"
    #Table 3 (unclassified sequences from stage 3) "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt"
    #Table 4 (no-blast-hits sequences from stage 3) "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt"
    concatenate_tables.pl "~{fasta_basename}_megablast_sample_specific_LCA.txt_classified.txt" "~{fasta_basename}_megablast_sample_specific_db_megablast_nt_LCA_classified.txt" "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_unclassified.txt" "~{fasta_basename}_megablast_sample_specific_db_megablast_nt.out_LCA.txt_no_hits.txt" > blastoff.txt
    awk 'BEGIN {FS="\t";OFS="\t"} {for (i=2; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' blastoff.txt > blastoff_no_first_column.txt
    mv blastoff_no_first_column.txt blastoff.txt

    LCA_table_to_kraken_output_format.pl blastoff.txt "~{trimmed_fasta}" > blastoff_krona.txt

    >>>
    output{

        File  most_popular_taxon_id = "sample_specific_db_taxa.txt"
        File  blastoff_results = "blastoff.txt"
        File  blastoff_krona = "blastoff_krona.txt"

    }
    runtime{
        docker:docker
        memory: machine_mem_gb + "GB"
        cpu: cpu
        disks: "local-disk" + disk_size_gb + "HDD"
        dx_instance_type: "n2-highmem-4"
    }

}
