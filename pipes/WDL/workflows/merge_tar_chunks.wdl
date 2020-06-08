version 1.0

import "../tasks/tasks_demux.wdl" as demux

workflow merge_tar_chunks {
    call demux.merge_tarballs
    output {
        File combined_tar = merge_tarballs.combined_tar
    }
}
