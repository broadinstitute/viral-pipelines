import "../tasks/tasks_read_utils.wdl" as read_utils

workflow merge_bams {
    call read_utils.merge_bams
}
