version 1.0

workflow Subsampler {
    input {
        File    metadata
        File    case_data
        File    keep_file
        File    remove_file
        File    filter_file

        String  id_column
        String  geo_column
        String  date_column
        String  baseline        =   "0.001"
        String  refgenome_size  =   "1"
        String  max_missing     =   "99"
        String  seed_num        =   "2007"
        String  start_date      =   "2020-03-01"
        String  end_date        =   "2021-12-31"
        String  unit            =   "week"

    }

    call subsample {
        input:
            metadata        =   metadata,
            case_data       =   case_data,
            keep_file       =   keep_file,
            remove_file     =   remove_file,
            filter_file     =   filter_file,
            id_column       =   id_column,
            geo_column      =   geo_column,
            date_column     =   date_column,
            baseline        =   baseline,
            refgenome_size  =   refgenome_size,
            max_missing     =   max_missing,
            seed_num        =   seed_num,
            start_date      =   start_date,
            end_date        =   end_date,
            unit            =   unit
    }

    output {
        File    matrix_days_file    =   subsample.genome_matrix_days
    }
}

task subsample {
    meta {
        description: "Run subsampler to get downsampled dataset and metadata."
    }
    input {
        File    metadata
        File    case_data
        File    keep_file
        File    remove_file
        File    filter_file

        String  id_column
        String  geo_column
        String  date_column
        String  baseline        =   "0.001"
        String  refgenome_size  =   "1"
        String  max_missing     =   "99"
        String  seed_num        =   "2007"
        String  start_date      =   "2020-03-01"
        String  end_date        =   "2021-12-31"
        String  unit            =   "week"

        String  docker = "quay.io/broadinstitute/subsampler"
    }
    command {

      # run the snakemake command
      # cd ../
      # /opt/conda/bin/snakemake subsample

      cd /opt/subsampler
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