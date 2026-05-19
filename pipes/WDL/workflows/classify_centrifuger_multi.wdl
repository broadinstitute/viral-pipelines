version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_centrifuger_multi {
    meta {
        description: "Batch Centrifuger taxonomic classification of multiple BAM files. All samples are classified on a single node with the database extracted once, amortizing the 200+ GB index load across all samples. Emits a per-read classification TSV and a Kraken2-style hierarchical kreport for each input BAM."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]+ reads_bams
        File         centrifuger_db_tgz
        String       db_name

        Int          machine_mem_gb = 240
        Int          cpu            = 8
    }

    parameter_meta {
        reads_bams: {
            description: "Array of BAM files, one per sample. Each BAM is converted to FASTQ and classified sequentially on a single node. Sample names are derived from BAM filenames via `basename <bam> .bam`.",
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
            description: "Memory in GB. Default 240 GB sized for multi-sample batch with NT-scale index."
        }
        cpu: {
            description: "Number of CPUs for centrifuger classify. Default 8."
        }
    }

    call metagenomics.centrifuger as run_centrifuger {
        input:
            reads_bams         = reads_bams,
            centrifuger_db_tgz = centrifuger_db_tgz,
            db_name            = db_name,
            machine_mem_gb     = machine_mem_gb,
            cpu                = cpu
    }

    output {
        Array[File] classification_tsvs = run_centrifuger.classification_tsvs
        Array[File] kreports            = run_centrifuger.kreports
        String      viralngs_version    = run_centrifuger.viralngs_version
    }
}
