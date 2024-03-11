version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_reports.wdl" as reports

workflow classify_krakenuniq {
    meta {
        description: "Taxonomic classification of reads using krakenuniq v1."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call metagenomics.krakenuniq

    call reports.aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = krakenuniq.krakenuniq_summary_reports
    }

    output {
        File        krakenuniq_krona_merged     = krakenuniq.krona_report_merged_html
        File        metagenomics_summary        = metag_summary_report.krakenuniq_aggregate_taxlevel_summary
        Array[File] krakenuniq_classified_reads = krakenuniq.krakenuniq_classified_reads
        Array[File] krakenuniq_summary_reports  = krakenuniq.krakenuniq_summary_reports
        Array[File] krakenuniq_krona_by_sample  = krakenuniq.krona_report_html
        String      viral_classify_version      = krakenuniq.viralngs_version
    }
}
