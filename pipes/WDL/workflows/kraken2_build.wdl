version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow kraken2_build {
    call metagenomics.build_kraken2_db

    output {
        File kraken2_db   = build_kraken2_db.kraken2_db
        File taxdump_tgz  = build_kraken2_db.taxdump_tgz
        File krona_db     = build_kraken2_db.krona_db
    }
}
