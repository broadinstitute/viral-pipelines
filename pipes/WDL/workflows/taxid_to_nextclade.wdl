version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow taxid_to_nextclade {
    meta {
        description: "Convert taxids to a nextclade dataset name"
    }

    call nextstrain.taxid_to_nextclade_dataset_name

    output {
        String nextclade_dataset    = taxid_to_nextclade_dataset_name.nextclade_dataset_name
    }
}
