version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow subsample_by_metadata {
    meta {
        description: "Filter and subsample a sequence set. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/filter.html"
    }

    call nextstrain.filter_subsample_sequences
    output {
        File filtered_fasta = filter_subsample_sequences.filtered_fasta
    }
}
