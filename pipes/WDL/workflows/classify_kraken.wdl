import "tasks_metagenomics.wdl" as metagenomics
import "tasks_reports.wdl" as reports

workflow classify_kraken {
    call metagenomics.kraken

    call reports.aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = kraken.kraken_summary_reports
    }
}
