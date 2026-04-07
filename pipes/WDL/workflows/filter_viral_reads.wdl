version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow filter_viral_reads {
    meta {
        description: "Filter a joined read-classification parquet to retain reads with any viral signal (K2 Viruses kingdom, Kallisto hit, VNP viral call, or geNomad viral taxonomy)."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File    joined_parquet
        String? sample_id_override
    }

    call metagenomics.filter_viral_reads as filter {
        input:
            joined_parquet     = joined_parquet,
            sample_id_override = sample_id_override
    }

    output {
        File viral_reads_parquet = filter.viral_reads_parquet
        File viral_read_ids      = filter.viral_read_ids
    }
}
