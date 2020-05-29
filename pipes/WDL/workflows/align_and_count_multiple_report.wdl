version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow align_and_count_multiple_report {
    meta {
        description: "Count the number of times reads map to provided reference sequences. Useful for counting spike-ins, etc."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+ reads_unmapped_bams
        File ref_db
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

    scatter(raw_reads in reads_unmapped_bams) {
        call reports.align_and_count {
            input:
                reads_bam = raw_reads,
                ref_db    = ref_db
        }
    }

    call reports.align_and_count_summary {
        input:
            counts_txt = align_and_count.report
    }

    call reports.align_and_count_summary as align_and_count_summary_top_hits {
        input:
            counts_txt    = align_and_count.report_top_hits,
            output_prefix = "count_summary_top_hits"
    }

    output {
        File report               = align_and_count_summary.count_summary
        File report_top_hits      = align_and_count_summary_top_hits.count_summary
        String viral_core_version = align_and_count_summary.viralngs_version
    }
}
