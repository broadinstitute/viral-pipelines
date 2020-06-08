version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_kaiju {
    call metagenomics.kaiju
    output {
        File    kaiju_report           = kaiju.kaiju_report
        File    kaiju_reads            = kaiju.kaiju_reads
        File    krona_report_html      = kaiju.krona_report_html
        String  viral_classify_version = kaiju.viralngs_version
    }
}
