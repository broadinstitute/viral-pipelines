version 1.0

import "../tasks/tasks_interhost.wdl" as interhost

workflow beast_gpu {
    call interhost.beast
    output {
        File        beast_log    = beast.beast_log
        Array[File] trees        = beast.trees
        String      beast_stdout = beast.beast_stdout
    }
}
