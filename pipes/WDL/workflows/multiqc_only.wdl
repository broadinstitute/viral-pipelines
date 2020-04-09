import "../tasks/tasks_reports.wdl" as reports

workflow multiqc_only {
    call reports.MultiQC
}
