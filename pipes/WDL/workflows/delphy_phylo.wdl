version 1.0

import "../tasks/tasks_delphy.wdl" as delphy

workflow delphy_phylo {
    meta {
        description: "Runs Delphy"
        author:      "Broad Viral Genomics"
        email:       "viral-ngs@broadinstitute.org"
    }

    call delphy.delphy

    output {
        File   delphy_trees            = delphy.delphy_trees
        File   delphy_log              = delphy.delphy_log
        File   delphy_for_web_ui       = delphy.delphy_for_web_ui
        File   delphy_output_for_beast = delphy.delphy_output_for_beast
        File   delphy_stdout           = delphy.delphy_stdout
        
        String delphy_version          = delphy.delphy_version
    }
}
