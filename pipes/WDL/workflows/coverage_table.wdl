version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow coverage_table {
    call reports.coverage_report
}
