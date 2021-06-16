version 1.0

import "../tasks/tasks_sarscov2.wdl" as sarscov2
import "../tasks/tasks_utils.wdl" as utils

workflow sarscov2_sequencing_reports {
    meta {
        description: "Produce per-state and per-collaborator weekly reports of SARS-CoV-2 surveillance data."
    }

    input {
        Array[File] assembly_stats_tsvs
        String?     max_date
    }

    call utils.today

    call utils.tsv_join {
        input:
            input_tsvs   = assembly_stats_tsvs,
            id_col       = 'sample',
            out_basename = 'assembly_stats-cumulative-~{max_date}'
    }

    call sarscov2.sequencing_report {
        input:
            assembly_stats_tsv = tsv_join.out_tsv,
            max_date           = select_first([max_date, today.date])
    }

    output {
        File        assembly_stats_cumulative_tsv = tsv_join.out_tsv
        Array[File] sequencing_reports_pdfs       = sequencing_report.reports
        Array[File] sequencing_reports_xlsxs      = sequencing_report.sheets
        File        sequencing_reports_zip        = sequencing_report.all_zip
        File        sequencing_report_tsv         = sequencing_report.all_tsv
    }
}
