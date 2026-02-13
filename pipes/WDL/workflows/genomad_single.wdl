version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow genomad_single {
    meta {
        description: "Runs genomad end-to-end classification on assembled contigs to identify viruses and plasmids. Accepts a FASTA assembly and genomad database, outputs classification summaries and annotated sequences."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File   assembly_fasta
        File   genomad_db_tgz

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
        cleanup: {
            description: "Delete intermediate genomad files after completion. Default true.",
            category: "common"
        }
    }

    call metagenomics.genomad_end_to_end {
        input:
            assembly_fasta = assembly_fasta,
            genomad_db_tgz = genomad_db_tgz,
            cleanup        = cleanup
    }

    call metagenomics.report_genomad_summary {
        input:
            virus_summary_tsv    = genomad_end_to_end.virus_summary,
            plasmid_summary_tsv  = genomad_end_to_end.plasmid_summary
    }

    # Extract optional counts/scores — truly undefined when zero
    if (length(report_genomad_summary.total_viruses_file) > 0) {
        Int extracted_total_viruses = read_int(report_genomad_summary.total_viruses_file[0])
    }
    if (length(report_genomad_summary.total_plasmids_file) > 0) {
        Int extracted_total_plasmids = read_int(report_genomad_summary.total_plasmids_file[0])
    }
    if (length(report_genomad_summary.top_virus_score_file) > 0) {
        Float extracted_top_virus_score = read_float(report_genomad_summary.top_virus_score_file[0])
    }

    output {
        # Genomad result files (always present, may be empty/header-only)
        File   genomad_virus_summary_tsv    = genomad_end_to_end.virus_summary
        File   genomad_plasmid_summary_tsv  = genomad_end_to_end.plasmid_summary
        File   genomad_virus_fasta          = genomad_end_to_end.virus_fasta
        File   genomad_plasmid_fasta        = genomad_end_to_end.plasmid_fasta

        # Summary statistics for Terra data tables
        Int?   genomad_total_viruses         = extracted_total_viruses
        Int?   genomad_total_plasmids        = extracted_total_plasmids
        Float? genomad_top_virus_score       = extracted_top_virus_score
        String genomad_top_virus_taxonomy    = report_genomad_summary.top_virus_taxonomy
    }
}
