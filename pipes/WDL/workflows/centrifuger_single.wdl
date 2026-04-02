version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow centrifuger_single {
    meta {
        description: "Centrifuger taxonomic classification of a single BAM file. Wraps the centrifuger task for single-sample use, suitable for Terra per-sample scatter."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File    reads_bam
        File    centrifuger_db_tgz
        String  db_name

        Int     machine_mem_gb = 240
        String  docker = "ghcr.io/broadinstitute/docker-centrifuger:1.0.0"
    }

    parameter_meta {
        reads_bam: {
            description: "Reads in BAM format for a single sample. Converted to FASTQ internally via picard SamToFastq.",
            patterns: ["*.bam"]
        }
        centrifuger_db_tgz: {
            description: "Pre-built Centrifuger index as a compressed tarball (.tar.gz, .tar.lz4, or .tar.zst).",
            patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.zst", "*.tar.bz2"]
        }
        db_name: {
            description: "Centrifuger index prefix (common filename stem of .1.cfr/.2.cfr/.3.cfr/.4.cfr files inside the tarball)."
        }
    }

    call metagenomics.centrifuger as run_centrifuger {
        input:
            reads_bams         = [reads_bam],
            centrifuger_db_tgz = centrifuger_db_tgz,
            db_name            = db_name,
            machine_mem_gb     = machine_mem_gb,
            docker             = docker
    }

    output {
        Array[File] classification_tsvs = run_centrifuger.classification_tsvs
        Array[File] kreports            = run_centrifuger.kreports
        Array[File] centrifuger_logs    = run_centrifuger.centrifuger_logs
    }
}
