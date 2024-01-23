version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow align_and_count_report {
    meta {
        description: "Align reads to reference with minimap2 and count the number of hits. Results are returned in the format of 'samtools idxstats'."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }
    parameter_meta {
        reads_unmapped_bams: {
            description: "Unaligned reads in BAM format",
            patterns: ["*.bam"]
        }
        ref_db: {
            description: "File containing sequences against which reads should me aligned and counted",
            patterns: ["*.fasta","*.fa"]
        }
    }

    call reports.align_and_count
    output {
        File   report             = align_and_count.report
        File   report_top_hits    = align_and_count.report_top_hits
        String viral_core_version = align_and_count.viralngs_version
    }
}
