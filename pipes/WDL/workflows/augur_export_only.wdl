version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow augur_export_only {
    meta {
        description: "Convert a newick formatted phylogenetic tree with other config settings and node values into a json suitable for auspice visualization. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/export.html"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call nextstrain.export_auspice_json
    output {
        File auspice_json       = export_auspice_json.virus_json
        File root_sequence_json = export_auspice_json.root_sequence_json
    }
}

