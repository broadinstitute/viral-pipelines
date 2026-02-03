version 1.1

import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_utils.wdl" as utils

workflow subsample_by_metadata_with_focal {
    meta {
        description: "Filter and subsample a global sequence set with a bias towards a geographic area of interest."
    }

    parameter_meta {
        sample_metadata_tsvs: {
            description: "Tab-separated metadata file that contain binning variables and values. Must contain all samples: output will be filtered to the IDs present in this file.",
            patterns: ["*.txt", "*.tsv"]
        }
        sequences_fasta: {
            description: "Sequences in fasta format.",
            patterns: ["*.fasta"]
        }

        focal_variable: {
            description: "The dataset will be bifurcated based on this column header."
        }
        focal_value: {
            description: "The dataset will be bifurcated based whether the focal_variable column matches this value or not. Rows that match this value are considered to be part of the 'focal' set of interest, rows that do not are part of the 'global' set."
        }

        focal_bin_variable: {
            description: "The focal subset of samples will be evenly subsampled across the discrete values of this column header."
        }
        focal_bin_max: {
            description: "The output will contain no more than this number of focal samples from each discrete value in the focal_bin_variable column."
        }

        global_bin_variable: {
            description: "The global subset of samples will be evenly subsampled across the discrete values of this column header."
        }
        global_bin_max: {
            description: "The output will contain no more than this number of global samples from each discrete value in the global_bin_variable column."
        }
    }

    input {
        Array[File]+ sample_metadata_tsvs
        File    sequences_fasta
        File?   priorities

        String  focal_variable      = "region"
        String  focal_value         = "North America"
        
        String  focal_bin_variable  = "division"
        Int     focal_bin_max       = 50
        
        String  global_bin_variable = "country"
        Int     global_bin_max      = 50
    }

    if(length(sample_metadata_tsvs)>1) {
        call utils.tsv_join {
            input:
                input_tsvs   = sample_metadata_tsvs,
                id_col       = 'strain',
                out_basename = "metadata-merged"
        }
    }

    call nextstrain.derived_cols {
        input:
            metadata_tsv = select_first(flatten([[tsv_join.out_tsv], sample_metadata_tsvs]))
    }

    call nextstrain.filter_subsample_sequences as prefilter {
        input:
            sequences_fasta     = sequences_fasta,
            sample_metadata_tsv = derived_cols.derived_metadata
    }

    call nextstrain.filter_subsample_sequences as subsample_focal {
        input:
            sequences_fasta     = prefilter.filtered_fasta,
            sample_metadata_tsv = derived_cols.derived_metadata,
            exclude_where       = ["${focal_variable}!=${focal_value}"],
            sequences_per_group = focal_bin_max,
            group_by            = focal_bin_variable,
            priority            = priorities
    }

    call nextstrain.filter_subsample_sequences as subsample_global {
        input:
            sequences_fasta     = prefilter.filtered_fasta,
            sample_metadata_tsv = derived_cols.derived_metadata,
            exclude_where       = ["${focal_variable}=${focal_value}"],
            sequences_per_group = global_bin_max,
            group_by            = global_bin_variable,
            priority            = priorities
    }

    call utils.concatenate as cat_fasta {
        input:
            infiles = [
                subsample_focal.filtered_fasta, subsample_global.filtered_fasta
            ],
            output_name = "subsampled.fasta"
    }

    call utils.fasta_to_ids {
        input:
            sequences_fasta = cat_fasta.combined
    }

    output {
        File metadata_merged      = derived_cols.derived_metadata
        File keep_list            = fasta_to_ids.ids_txt
        File subsampled_sequences = cat_fasta.combined
        Int  focal_kept           = subsample_focal.sequences_out
        Int  global_kept          = subsample_global.sequences_out
        Int  sequences_kept       = subsample_focal.sequences_out + subsample_global.sequences_out
    }
}
