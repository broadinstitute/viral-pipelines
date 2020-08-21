version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow augur_from_beast_mcc {
    meta {
        description: "Visualize BEAST output with Nextstrain. This workflow converts a BEAST MCC tree (.tree file) into an Auspice v2 json file. See https://nextstrain-augur.readthedocs.io/en/stable/faq/import-beast.html for details."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        File    beast_mcc_tree
        File    auspice_config
    }

    parameter_meta {
        beast_mcc_tree: {
          description: "A maximum clade credibility (MCC) tree (.tree file) that is output from a BEAST run.",
          patterns: ["*.tree"]
        }
        auspice_config: {
          description: "A file specifying options to customize the auspice export; see: https://nextstrain.github.io/auspice/customise-client/introduction",
          patterns: ["*.json", "*.txt"]
        }
    }

    call nextstrain.augur_import_beast {
        input:
            beast_mcc_tree = beast_mcc_tree
    }
    call nextstrain.export_auspice_json {
        input:
            tree            = augur_import_beast.tree_newick,
            node_data_jsons = [augur_import_beast.node_data_json],
            auspice_config  = auspice_config
    }

    output {
        File  beast_mcc_tree_newick      = augur_import_beast.tree_newick
        File  auspice_input_json         = export_auspice_json.virus_json
    }
}