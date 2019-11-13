import "tasks_demux.wdl" as demux

workflow merge_tarballs {
    call demux.merge_tarballs
}
