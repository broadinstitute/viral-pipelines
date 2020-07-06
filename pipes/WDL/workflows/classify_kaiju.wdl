version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_kaiju {
    meta {
        description: "Taxonomic classification of reads with kaiju."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    call metagenomics.kaiju
    output {
        File    kaiju_report           = kaiju.kaiju_report
        File    kaiju_reads            = kaiju.kaiju_reads
        File    krona_report_html      = kaiju.krona_report_html
        String  viral_classify_version = kaiju.viralngs_version
    }
}
