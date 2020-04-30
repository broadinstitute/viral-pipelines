version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow spikein {
    call reports.spikein_report
    output {
        File report               = spikein_report.report
        File report_top_hits      = spikein_report.report_top_hits
        String viral_core_version = spikein_report.viralngs_version
    }
}
