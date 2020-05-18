version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow newick_to_auspice {
    meta {
        description: "Convert a newick formatted phylogenetic tree into a json suitable for auspice visualization. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/export.html"
    }

    call nextstrain.export_auspice_json
    output {
        File auspice_json = export_auspice_json.virus_json
    }
}
