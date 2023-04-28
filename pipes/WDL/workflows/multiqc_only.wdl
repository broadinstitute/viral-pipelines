version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow multiqc_only {
    meta {
        description: "Combine multiple FastQC reports into a single MultiQC summary."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }
    input {
        Array[File] input_files
        String      file_name = "multiqc-raw.html"
    }
    call reports.MultiQC {
        input:
            input_files = input_files,
            file_name = file_name
    }
    output {
        File multiqc = MultiQC.multiqc_report
    }
}
