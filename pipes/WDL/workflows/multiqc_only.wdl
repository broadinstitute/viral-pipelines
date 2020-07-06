version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow multiqc_only {
    meta {
        description: "Combine multiple FastQC reports into a single MultiQC summary."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }
    call reports.MultiQC
    output {
        File        multiqc = MultiQC.multiqc_report
    }
}
