version 1.0

import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_read_utils.wdl" as read_utils

workflow align_and_generate_PAF {
    meta {
        description: "Aligns and generates a report of alignment results (in PAF format) for a set of unmapped BAM files against a reference database."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    call assembly.align_reads as align

    call read_utils.BamToPAF {
        input:
            aligned_bam = align.aligned_bam
    }

    output {
        File   aligned_bam                   = align.aligned_bam
        File   aligned_bam_idx               = align.aligned_bam_idx
        File   aligned_bam_flagstat          = align.aligned_bam_flagstat

        File   aligned_PAF                   = BamToPAF.aligned_paf
        
        String viral_core_version            = align.viralngs_version
    }
}