import "../tasks/tasks_reports.wdl" as reports

workflow fastqc_only {
    call reports.fastqc
}
