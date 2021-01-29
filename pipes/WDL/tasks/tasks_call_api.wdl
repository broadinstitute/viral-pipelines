version 1.0

task  upload_entities_tsv{
  input {
    String        workspace_name
    String        terra_project
    File          tsv_file

    String        docker="schaluvadi/pathogen-genomic-surveillance:api-wdl"
  }

  command {

    python3 /projects/cdc-sabeti-covid-19/create_data_tables.py -t ~{tsv_file} \
            -p ~{terra_project} \
            -w ~{workspace_name}

  }

  runtime {
    docker: docker
    preemptible: 0
  }

  output {
    String  status  = read_string(stdout())
  }
}