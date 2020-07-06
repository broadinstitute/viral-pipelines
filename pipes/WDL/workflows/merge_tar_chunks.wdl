version 1.0

import "../tasks/tasks_demux.wdl" as demux

workflow merge_tar_chunks {
    meta {
        description: "Combine multiple tar files (possibly compressed by gzip, bz2, lz4, zstd, etc) into a single tar file. Originally meant for combining streaming upload chunks from a sequencing run."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call demux.merge_tarballs
    output {
        File combined_tar = merge_tarballs.combined_tar
    }
}
