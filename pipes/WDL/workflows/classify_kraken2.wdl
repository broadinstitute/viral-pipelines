version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_kraken2 {
    call metagenomics.kraken2

    output {
        File    kraken2_reads_report   = kraken2.kraken2_reads_report
        File    kraken2_summary_report = kraken2.kraken2_summary_report
        File    krona_report_html      = kraken2.krona_report_html
        String  viral_classify_version = kraken2.viralngs_version
    }
}
