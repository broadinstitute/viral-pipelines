version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_reports.wdl" as reports

workflow merge_metagenomics {
    meta {
        description: "Combine metagenomic reports from single samples into an aggregate report."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File] krakenuniq_summary_reports
        Array[File] krona_per_sample
    }

    call metagenomics.krona_merge {
        input:
            krona_reports = krona_per_sample,
            out_basename = "merged_krona"
    }

    call reports.aggregate_metagenomics_reports {
        input:
            kraken_summary_reports = krakenuniq_summary_reports
    }

    output {
        File krona_merged         = krona_merge.krona_report_html
        File metagenomics_summary = aggregate_metagenomics_reports.krakenuniq_aggregate_taxlevel_summary
    }
}
