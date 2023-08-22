version 1.0

workflow Subsampler {
    input {
        File    metadata
        String  geo_column
        String  date_column
    }

    call subsample {
        input:
            metadata    =   metadata,
            geo_column  =   geo_column,
            date_column =   date_column
    }
}

task subsample {
    meta {
        description: "Run subsampler to get downsampled dataset and metadata."
    }
    input {
        File    metadata
        String  geo_column
        String  date_column
        String  docker = "quay.io/broadinstitute/subsampler"
    }
    command {

      # run the snakemake command
      # cd ../
      # /opt/conda/bin/snakemake subsample

      snakemake genome_matrix --config metadata=~{metadata} \
                                       geo_column=~{geo_column} \
                                       date_column=~{date_column}

    }
    runtime {
        docker: docker
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File genome_matrix_days = "genome_matrix_days.tsv"
    }
}