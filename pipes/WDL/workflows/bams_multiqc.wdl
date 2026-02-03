version 1.1

import "../tasks/tasks_reports.wdl" as reports

workflow bams_multiqc {
    meta {
        description: "Run FastQC on a set of BAM files, and then MultiQC to summarize all outputs."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]+ read_bams
    }

    call reports.multiqc_from_bams {
        input:
            input_bams = read_bams
    }

    output {
        File        multiqc  = multiqc_from_bams.multiqc_report
        Array[File] fastqcs  = multiqc_from_bams.fastqc_html
    }
}
