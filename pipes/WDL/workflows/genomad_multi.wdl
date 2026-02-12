version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow genomad_multi {
    meta {
        description: "Runs genomad end-to-end classification on multiple assemblies in parallel, identifying viruses and plasmids in each sample."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]+ assembly_fastas
        File         genomad_db_tgz

        Boolean      cleanup = true
    }

    parameter_meta {
        assembly_fastas: {
            description: "Array of assembled contigs in FASTA format to classify for viral and plasmid sequences.",
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

    scatter(assembly_fasta in assembly_fastas) {
        call metagenomics.genomad_end_to_end {
            input:
                assembly_fasta = assembly_fasta,
                genomad_db_tgz = genomad_db_tgz,
                cleanup        = cleanup
        }
    }

    output {
        Array[File] virus_summary_tsvs   = genomad_end_to_end.virus_summary
        Array[File] plasmid_summary_tsvs = genomad_end_to_end.plasmid_summary
        Array[File] virus_fastas         = genomad_end_to_end.virus_fasta
        Array[File] plasmid_fastas       = genomad_end_to_end.plasmid_fasta
        Array[Int]  genomad_max_ram_gb   = genomad_end_to_end.max_ram_gb
        String      viral_classify_version = genomad_end_to_end.viralngs_version[0]
    }
}
