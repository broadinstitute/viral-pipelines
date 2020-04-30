version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow align_and_count_report {
    call reports.align_and_count
    output {
        File report               = align_and_count.report
        File report_top_hits      = align_and_count.report_top_hits
        String viral_core_version = align_and_count.viralngs_version
    }
}
