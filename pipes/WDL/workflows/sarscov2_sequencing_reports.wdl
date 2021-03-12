version 1.0

import "../tasks/tasks_sarscov2.wdl" as sarscov2

workflow sarscov2_sequencing_reports {
    meta {
        description: "Produce per-state and per-collaborator weekly reports of SARS-CoV-2 surveillance data."
    }

    call sarscov2.sequencing_report 

    output {
        Array[File]  sequencing_reports_pdfs  = sequencing_report.reports
        Array[File]  sequencing_reports_xlsxs = sequencing_report.sheets
    }
}