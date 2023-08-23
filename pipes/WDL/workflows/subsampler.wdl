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
        String? start_date
        String? end_date
        String  baseline        =   "0.001"
        String  refgenome_size  =   "1"
        String  max_missing     =   "99"
        String  seed_num        =   "2007"
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
        File    genome_matrix_days_file             =   subsample.genome_matrix_days
        File    matrix_genomes_unit_file            =   subsample.matrix_genomes_unit
        File    matrix_cases_unit_file              =   subsample.matrix_cases_unit
        File    weekly_sampling_proportions_file    =   subsample.weekly_sampling_proportions
		File    weekly_sampling_bias_file           =   subsample.weekly_sampling_bias
		File    matrix_genomes_unit_corrected_file  =   subsample.matrix_genomes_unit_corrected
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

        String? start_date
        String? end_date
        String  id_column
        String  geo_column
        String  date_column
        String  baseline        =   "0.001"
        String  refgenome_size  =   "1"
        String  max_missing     =   "99"
        String  seed_num        =   "2007"
        String  unit            =   "week"

        String  docker = "quay.io/broadinstitute/subsampler"
    }
    command {
        # run the snakemake command
        cd /opt/subsampler
        
        # snakemake genome_matrix --config metadata=~{metadata} \
        #                                  case_data=~{case_data} \
        #                                  keep_file=~{keep_file} \
        #                                  remove_file=~{remove_file} \
        #                                  filter_file=~{filter_file} \
        #                                  id_column=~{id_column} \
        #                                  geo_column=~{geo_column} \
        #                                  date_column=~{date_column} \
        #                                  baseline=~{baseline} \
        #                                  refgenome_size=~{refgenome_size} \
        #                                  max_missing=~{max_missing} \
        #                                  seed_num=~{seed_num} \
        #                                  start_date=~{start_date} \
        #                                  end_date=~{end_date} \
        #                                  unit=~{unit}

        # snakemake unit_conversion --config metadata=~{metadata} \
        #                                  case_data=~{case_data} \
        #                                  keep_file=~{keep_file} \
        #                                  remove_file=~{remove_file} \
        #                                  filter_file=~{filter_file} \
        #                                  id_column=~{id_column} \
        #                                  geo_column=~{geo_column} \
        #                                  date_column=~{date_column} \
        #                                  baseline=~{baseline} \
        #                                  refgenome_size=~{refgenome_size} \
        #                                  max_missing=~{max_missing} \
        #                                  seed_num=~{seed_num} \
        #                                  start_date=~{start_date} \
        #                                  end_date=~{end_date} \
        #                                  unit=~{unit}

        # snakemake correct_bias --config metadata=~{metadata} \
        #                                  case_data=~{case_data} \
        #                                  keep_file=~{keep_file} \
        #                                  remove_file=~{remove_file} \
        #                                  filter_file=~{filter_file} \
        #                                  id_column=~{id_column} \
        #                                  geo_column=~{geo_column} \
        #                                  date_column=~{date_column} \
        #                                  baseline=~{baseline} \
        #                                  refgenome_size=~{refgenome_size} \
        #                                  max_missing=~{max_missing} \
        #                                  seed_num=~{seed_num} \
        #                                  start_date=~{start_date} \
        #                                  end_date=~{end_date} \
        #                                  unit=~{unit}

        snakemake subsample --config metadata=~{metadata} \
                                         case_data=~{case_data} \
                                         keep_file=~{keep_file} \
                                         remove_file=~{remove_file} \
                                         filter_file=~{filter_file} \
                                         id_column=~{id_column} \
                                         geo_column=~{geo_column} \
                                         date_column=~{date_column} \
                                         baseline=~{baseline} \
                                         refgenome_size=~{refgenome_size} \
                                         max_missing=~{max_missing} \
                                         seed_num=~{seed_num} \
                                         unit=~{unit} \
                                         ~{"start_date=" + start_date} \
                                         ~{"end_date=" + end_date}
        # delocalization script requires outputs to be in /cromwell_root
        cp outputs/* /cromwell_root
        # cp *.txt /cromwell_root
        cd /cromwell_root
    }
    runtime {
        docker: docker
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File    genome_matrix_days              =   "genome_matrix_days.tsv"
        File    matrix_genomes_unit             =   "matrix_genomes_unit.tsv"
        File    matrix_cases_unit               =   "matrix_cases_unit.tsv"
        File    weekly_sampling_proportions     =   "weekly_sampling_proportions.tsv"
		File    weekly_sampling_bias            =   "weekly_sampling_bias.tsv"
		File    matrix_genomes_unit_corrected   =   "matrix_genomes_unit_corrected.tsv"
    }
}