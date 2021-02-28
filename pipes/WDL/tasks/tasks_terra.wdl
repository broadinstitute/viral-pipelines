version 1.0


task upload_entities_tsv {
  input {
    String        workspace_name
    String        terra_project
    File          tsv_file
    Array[String] cleaned_reads_unaligned_bams_string
    File          meta_by_filename_json

    String        docker="schaluvadi/pathogen-genomic-surveillance:api-wdl"
    Int?          machine_mem_gb
  }

  command {
    
    echo ~{sep="," cleaned_reads_unaligned_bams_string} > cleaned_bam_strings.txt

    python3 /projects/cdc-sabeti-covid-19/create_data_tables.py -t ~{tsv_file} \
            -p ~{terra_project} \
            -w ~{workspace_name} \
            -b cleaned_bam_strings.txt \
            -j ~{meta_by_filename_json}

  }

  runtime {
    docker: docker
    preemptible: 0
    memory: select_first([machine_mem_gb, 2]) + " GB"
  }

  output {
    String  status  = read_string(stdout())
  }
}