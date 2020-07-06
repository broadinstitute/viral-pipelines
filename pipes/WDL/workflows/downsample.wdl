version 1.0

import "../tasks/tasks_read_utils.wdl" as reads

workflow downsample {
    meta {
        description: "Random subsampling of reads."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call reads.downsample_bams

    output {
      Array[File] downsampled_bam  = downsample_bams.downsampled_bam
      String      viralngs_version = downsample_bams.viralngs_version
    }
}
