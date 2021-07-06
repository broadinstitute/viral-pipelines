version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow sarscov2_gisaid_ingest {
    meta {
        description: "Sanitize data downloaded from GISAID for use in Nextstrain/augur. See: https://nextstrain.github.io/ncov/data-prep#curate-data-from-the-full-gisaid-database"
    }
    input {
        File sequences_gisaid_fasta
        File metadata_gisaid_tsv

        String? gcs_out
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
        gcs_out: {
            description: "If specified, GCP bucket prefix for storage of the output data."
        }
    }
    call nextstrain.nextstrain_ncov_sanitize_gisaid_data as sanitize_gisaid {
        input:
            sequences_gisaid_fasta = sequences_gisaid_fasta,
            metadata_gisaid_tsv    = metadata_gisaid_tsv
    }

    if(defined(gcs_out)) {
        call terra.gcs_copy as gcs_dump {
            input:
                infiles        = [sanitize_gisaid.metadata_gisaid_sanitized_tsv, sanitize_gisaid.sequences_gisaid_sanitized_fasta],
                gcs_uri_prefix = "~{gcs_out}/"
        }
    }

    output {
        File sequences_gisaid_sanitized_fasta = sanitize_gisaid.sequences_gisaid_sanitized_fasta
        File metadata_gisaid_sanitized_tsv    = sanitize_gisaid.metadata_gisaid_sanitized_tsv
    }
}

