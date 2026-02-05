version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow coverage_table {
    call reports.coverage_report
    output {
        File   coverage_report_txt = coverage_report.coverage_report
        String viral_core_version  = coverage_report.viralngs_version
    }
}
