version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow sarscov2_gisaid_ingest {
    meta {
        description: "Sanitize data downloaded from GISAID for use in Nextstrain/augur. See: https://nextstrain.github.io/ncov/data-prep#curate-data-from-the-full-gisaid-database"
    }
    input {
        File sequences_gisaid_fasta
        File metadata_gisaid_tsv
    }
    parameter_meta {
        sequences_gisaid_fasta: {
          description: "Multiple sequences downloaded from GISAID",
          patterns: ["*.fasta","*.fasta.xy","*.fasta.gz"]
        }
        metadata_gisaid_tsv: {
            description: "Tab-separated metadata file for sequences downloaded from GISAID and passed in via sequences_gisaid_fasta.",
            patterns: ["*.txt", "*.tsv","*.tsv.xy","*.tsv.gz"]
        }
    }
    call nextstrain.nextstrain_ncov_sanitize_gisaid_data as sanitize_gisaid {
        input:
            sequences_gisaid_fasta = sequences_gisaid_fasta,
            metadata_gisaid_tsv    = metadata_gisaid_tsv
    }
    output {
        File sequences_gisaid_sanitized_fasta    = sanitize_gisaid.sequences_gisaid_sanitized_fasta
        File metadata_gisaid_sanitized_tsv = sanitize_gisaid.metadata_gisaid_sanitized_tsv
    }
}

