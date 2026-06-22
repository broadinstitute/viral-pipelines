version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_genomad {
    meta {
        description: "Runs geNomad end-to-end classification on assembled contigs to identify viruses, plasmids, and proviruses."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File   assembly_fasta
        File   genomad_db_tgz

        Boolean cleanup = true
    }

    parameter_meta {
        assembly_fasta: {
            description: "Assembled contigs in FASTA format to classify for viral, plasmid, and proviral sequences.",
            patterns: ["*.fasta", "*.fa", "*.fna"],
            category: "required"
        }
        genomad_db_tgz: {
            description: "Pre-built geNomad database tarball. The committed test input uses the public Zenodo URL so clean checkouts do not require the ignored local test/input/genomad_db-test.tar.zst fixture; local runs may override this with a local database tarball.",
            patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"],
            category: "required"
        }
        cleanup: {
            description: "Delete intermediate geNomad files after completion. Default true.",
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
            plasmid_summary_tsv  = genomad_end_to_end.plasmid_summary,
            provirus_summary_tsv = genomad_end_to_end.provirus_summary
    }

    if (length(report_genomad_summary.total_viruses_file) > 0) {
        Int extracted_total_viruses = read_int(report_genomad_summary.total_viruses_file[0])
    }
    if (length(report_genomad_summary.total_plasmids_file) > 0) {
        Int extracted_total_plasmids = read_int(report_genomad_summary.total_plasmids_file[0])
    }
    if (length(report_genomad_summary.total_proviruses_file) > 0) {
        Int extracted_total_proviruses = read_int(report_genomad_summary.total_proviruses_file[0])
    }
    if (length(report_genomad_summary.top_virus_score_file) > 0) {
        Float extracted_top_virus_score = read_float(report_genomad_summary.top_virus_score_file[0])
    }

    output {
        File   genomad_virus_summary_tsv     = genomad_end_to_end.virus_summary
        File   genomad_plasmid_summary_tsv   = genomad_end_to_end.plasmid_summary
        File   genomad_provirus_summary_tsv  = genomad_end_to_end.provirus_summary
        File   genomad_virus_fasta           = genomad_end_to_end.virus_fasta
        File   genomad_plasmid_fasta         = genomad_end_to_end.plasmid_fasta
        File   genomad_proteins_faa          = genomad_end_to_end.proteins_faa

        Int    genomad_max_ram_gb            = genomad_end_to_end.max_ram_gb
        String viral_classify_version        = genomad_end_to_end.viralngs_version

        Int?   genomad_total_viruses         = extracted_total_viruses
        Int?   genomad_total_plasmids        = extracted_total_plasmids
        Int?   genomad_total_proviruses      = extracted_total_proviruses
        Float? genomad_top_virus_score       = extracted_top_virus_score
        String genomad_top_virus_taxonomy    = report_genomad_summary.top_virus_taxonomy
    }
}
