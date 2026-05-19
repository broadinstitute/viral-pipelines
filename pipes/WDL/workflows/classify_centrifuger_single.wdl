version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_centrifuger_single {
    meta {
        description: "Centrifuger taxonomic classification of a single BAM file. Wraps the centrifuger task for single-sample use, suitable for Terra per-sample scatter. Emits a per-read classification TSV and a Kraken2-style hierarchical kreport."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File         reads_bam
        File         centrifuger_db_tgz
        String       db_name

        Int          machine_mem_gb = 240
        Int          cpu            = 8
    }

    parameter_meta {
        reads_bam: {
            description: "Single unaligned BAM file. Sample name is derived from the filename via `basename <bam> .bam`.",
            patterns: ["*.bam"]
        }
        centrifuger_db_tgz: {
            description: "Pre-built Centrifuger index as a compressed tarball.",
            patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.zst", "*.tar.bz2"]
        }
        db_name: {
            description: "Centrifuger index prefix (common filename stem of .1.cfr/.2.cfr/.3.cfr/.4.cfr files inside the tarball)."
        }
        machine_mem_gb: {
            description: "Memory in GB. Default 240 GB sized for NT-scale centrifuger index."
        }
        cpu: {
            description: "Number of CPUs for centrifuger classify. Default 8."
        }
    }

    call metagenomics.centrifuger as run_centrifuger {
        input:
            reads_bams         = [reads_bam],
            centrifuger_db_tgz = centrifuger_db_tgz,
            db_name            = db_name,
            machine_mem_gb     = machine_mem_gb,
            cpu                = cpu
    }

    output {
        File   classification_tsv = run_centrifuger.classification_tsvs[0]
        File   kreport            = run_centrifuger.kreports[0]
        String viralngs_version   = run_centrifuger.viralngs_version
    }
}
