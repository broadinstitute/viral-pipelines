version 1.0

import "../tasks/tasks_interhost.wdl" as interhost

workflow beast_gpu {
    meta {
        description: "Runs BEAST (v1) on a GPU instance. Use with care--this can be expensive if run incorrectly."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call interhost.beast
    output {
        File        beast_log    = beast.beast_log
        Array[File] trees        = beast.trees
        String      beast_stdout = beast.beast_stdout
    }
}
