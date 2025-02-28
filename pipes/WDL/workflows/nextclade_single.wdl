version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow nextclade_single {
    meta {
        description: "Run Nextclade on a single genome"
    }

    call nextstrain.taxid_to_nextclade_dataset_name

    call nextstrain.nextclade_one_sample {
        input:
            dataset_name = taxid_to_nextclade_dataset_name.nextclade_dataset_name
    }

    output {
        String nextclade_clade    = nextclade_one_sample.nextclade_clade
        File   nextclade_tsv      = nextclade_one_sample.nextclade_tsv
        File   nextclade_json     = nextclade_one_sample.nextclade_json
        String nextclade_aa_subs  = nextclade_one_sample.aa_subs_csv
        String nextclade_aa_dels  = nextclade_one_sample.aa_dels_csv
        String nextclade_shortclade = nextclade_one_sample.nextclade_shortclade
        String nextclade_subclade = nextclade_one_sample.nextclade_subclade
        String nextclade_pango    = nextclade_one_sample.nextclade_pango
        String nextclade_version  = nextclade_one_sample.nextclade_version
    }
}
