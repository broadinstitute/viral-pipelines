version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow genomad_single {
    meta {
        description: "Runs genomad end-to-end classification on assembled contigs to identify viruses, plasmids, and proviruses. Accepts a FASTA assembly and genomad database, outputs classification summaries and annotated sequences."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File   assembly_fasta
        File   genomad_db_tgz

        Int    min_length = 2500
        Boolean cleanup   = true
    }

    parameter_meta {
        assembly_fasta: {
            description: "Assembled contigs in FASTA format to classify for viral and plasmid sequences.",
            patterns: ["*.fasta", "*.fa", "*.fna"],
            category: "required"
        }
        genomad_db_tgz: {
            description: "Pre-built genomad database tarball. Download from https://zenodo.org/records/10594875",
            patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"],
            category: "required"
        }
        min_length: {
            description: "Minimum sequence length (bp) for genomad classification. Default 2500.",
            category: "common"
        }
        cleanup: {
            description: "Delete intermediate genomad files after completion. Default true.",
            category: "common"
        }
    }

    call metagenomics.genomad_end_to_end {
        input:
            assembly_fasta = assembly_fasta,
            genomad_db_tgz = genomad_db_tgz,
            min_length     = min_length,
            cleanup        = cleanup
    }

    call metagenomics.report_genomad_summary {
        input:
            virus_summary_tsv    = genomad_end_to_end.virus_summary,
            plasmid_summary_tsv  = genomad_end_to_end.plasmid_summary,
            provirus_summary_tsv = genomad_end_to_end.provirus_summary
    }

    output {
        # Genomad result files (optional — absent when nothing found)
        File?  virus_summary_tsv    = genomad_end_to_end.virus_summary
        File?  plasmid_summary_tsv  = genomad_end_to_end.plasmid_summary
        File?  provirus_summary_tsv = genomad_end_to_end.provirus_summary
        File?  virus_fasta          = genomad_end_to_end.virus_fasta
        File?  plasmid_fasta        = genomad_end_to_end.plasmid_fasta

        # Summary statistics for Terra data tables
        Int    total_viruses         = report_genomad_summary.total_viruses
        Int    total_plasmids        = report_genomad_summary.total_plasmids
        Int    total_proviruses      = report_genomad_summary.total_proviruses
        String top_virus_name        = report_genomad_summary.top_virus_name
        Float  top_virus_score       = report_genomad_summary.top_virus_score
        String top_virus_taxonomy    = report_genomad_summary.top_virus_taxonomy

        # Runtime metadata
        Int    genomad_max_ram_gb    = genomad_end_to_end.max_ram_gb
        String viral_classify_version = genomad_end_to_end.viralngs_version
    }
}
