version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow parse_kraken2_reads {
    meta {
        description: "Parse Kraken2 per-read output and annotate each read with NCBI taxonomy (name, kingdom, rank) using a pre-built DuckDB taxonomy database."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File    kraken2_reads_output
        File    taxonomy_db
        String  sample_id_suffix_re = "\\.l0\\d+.*$"
        Boolean resolve_strains = false
    }

    call metagenomics.parse_kraken2_reads as parse_reads {
        input:
            kraken2_reads_output = kraken2_reads_output,
            taxonomy_db          = taxonomy_db,
            sample_id_suffix_re  = sample_id_suffix_re,
            resolve_strains      = resolve_strains
    }

    output {
        File read_taxonomy = parse_reads.read_taxonomy
    }
}
