import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_reports.wdl" as reports

workflow classify_krakenuniq {
    call metagenomics.krakenuniq

    call metagenomics.krona_merge {
        input:
            krona_reports = krakenuniq.krona_report_html,
            out_basename  = "krakenuniq.krona.combined.html"
    }

    call reports.aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = krakenuniq.krakenuniq_summary_reports
    }
}
