version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow beast_to_auspice {
    meta {
        description: "Visualize BEAST output with Nextstrain. This workflow converts a BEAST MCC tree (.tree file) into an Auspice v2 json file. See https://nextstrain-augur.readthedocs.io/en/stable/faq/import-beast.html for details."
    }

    input {
        File            beast_mcc_tree
        File            sample_metadata
        String          virus
    }

    parameter_meta {
        beast_mcc_tree: {
          description: "A maximum clade credibility (MCC) tree (.tree file) that is output from a BEAST run.",
          patterns: ["*.tree"]
        }
        sample_metadata: {
          description: "Metadata in tab-separated text format. See https://nextstrain-augur.readthedocs.io/en/stable/faq/metadata.html for details.",
          patterns: ["*.txt", "*.tsv"]
        }
        virus: {
          description: "A filename-friendly string that is used as a base for output file names."
        }
    }

    call nextstrain.augur_import_beast {
        input:
            beast_mcc_tree = beast_mcc_tree
    }
    call nextstrain.export_auspice_json {
        input:
            tree            = augur_import_beast.tree_newick,
            metadata        = sample_metadata,
            node_data_jsons = [augur_import_beast.node_data_json]
    }

    output {
        File  beast_mcc_tree_newick      = augur_import_beast.tree_newick
        File  node_data_json             = augur_import_beast.node_data_json
        File  auspice_input_json         = export_auspice_json.virus_json
    }
}
