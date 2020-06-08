version 1.0

import "../tasks/tasks_intrahost.wdl" as intrahost

workflow isnvs_one_sample {
    call intrahost.isnvs_per_sample

    output {
        File   isnvsFile                   = isnvs_per_sample.isnvsFile
        String isnvs_viral_phylo_version   = isnvs_per_sample.viralngs_version
    }
}
