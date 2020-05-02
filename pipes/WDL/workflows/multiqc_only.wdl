version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow multiqc_only {
    call reports.MultiQC
    output {
        File        multiqc = MultiQC.multiqc_report
    }
}
