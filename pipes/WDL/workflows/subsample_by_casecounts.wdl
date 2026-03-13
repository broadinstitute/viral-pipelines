version 1.1

import "../tasks/tasks_interhost.wdl" as interhost

workflow subsampler_only {

    call interhost.subsample_by_cases

    output {
        File    genome_matrix_days_file             =   subsample_by_cases.genome_matrix_days
        File    matrix_genomes_unit_file            =   subsample_by_cases.matrix_genomes_unit
        File    matrix_cases_unit_file              =   subsample_by_cases.matrix_cases_unit
        File    weekly_sampling_proportions_file    =   subsample_by_cases.weekly_sampling_proportions
		File    weekly_sampling_bias_file           =   subsample_by_cases.weekly_sampling_bias
		File    matrix_genomes_unit_corrected_file  =   subsample_by_cases.matrix_genomes_unit_corrected
        File    selected_sequences_file             =   subsample_by_cases.selected_sequences
        File    selected_metadata_file              =   subsample_by_cases.selected_metadata
        File    sampling_stats_file                 =   subsample_by_cases.sampling_stats
    }
}
