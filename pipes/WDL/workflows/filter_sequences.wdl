version 1.1

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow filter_sequences {
    meta {
        description: "Filter and subsample a sequence set. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/filter.html"
    }

    call nextstrain.filter_subsample_sequences as filter
}
