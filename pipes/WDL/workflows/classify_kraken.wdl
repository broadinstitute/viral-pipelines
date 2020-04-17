import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_reports.wdl" as reports

workflow classify_kraken {
    call metagenomics.kraken

    call metagenomics.krona_merge {
        input:
            krona_reports = kraken.krona_report_html,
            out_basename  = "kraken.krona.combined.html"
    }

    call reports.aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = kraken.kraken_summary_reports
    }
}
