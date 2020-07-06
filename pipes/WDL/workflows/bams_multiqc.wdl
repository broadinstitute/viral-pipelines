version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow bams_multiqc {
    meta {
        description: "Run FastQC on a set of BAM files, and then MultiQC to summarize all outputs."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+  read_bams
    }

    scatter(reads_bam in read_bams) {
        call reports.fastqc as fastqc {
            input:
                reads_bam = reads_bam
        }
    }

    call reports.MultiQC {
        input:
            input_files = fastqc.fastqc_zip
    }

    output {
        File        multiqc = MultiQC.multiqc_report
        Array[File] fastqcs = fastqc.fastqc_html
        String      viral_core_version = fastqc.viralngs_version[0]
    }
}
