version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow align_and_generate_reads_report {
    meta {
        description: "Aligns and generates a report of alignment results for a set of unmapped BAM files against a reference database. The report includes per-read mapping level statistics."
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

    call reports.align_and_generate_reads_report

    output {
        File   aligned_reads_paf     = align_and_generate_reads_report.aligned_paf
        File   aligned_reads_bam     = align_and_generate_reads_report.aligned_bam
        
        String viral_core_version    = align_and_generate_reads_report.viralngs_version
    }
}
