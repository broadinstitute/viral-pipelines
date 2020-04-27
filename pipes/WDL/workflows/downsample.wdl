version 1.0

import "../tasks/tasks_read_utils.wdl" as reads

workflow downsample {
    call reads.downsample_bams
}
