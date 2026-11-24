version 1.1

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow join_read_classifications {
    meta {
        description: "Annotates native Kraken2 and Centrifuger per-read output with a prebuilt taxonomy DuckDB, then joins Kallisto, Kraken2, Centrifuger, VirNucPro, and geNomad read-level calls into one Parquet row per sample and read."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    call metagenomics.join_read_classifications as join_reads

    output {
        File   read_classifications_parquet = join_reads.read_classifications_parquet
        Int    joined_read_count            = join_reads.joined_read_count
        String viralngs_version             = join_reads.viralngs_version
    }
}
