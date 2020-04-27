version 1.0

import "../tasks/tasks_read_utils.wdl" as read_utils

workflow merge_bams {
    call read_utils.merge_and_reheader_bams
}
