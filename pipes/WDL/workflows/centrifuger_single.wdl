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

        File?   taxonomy_db
        Boolean resolve_strains = false

        Int     machine_mem_gb = 256
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
        taxonomy_db: {
            description: "Pre-built DuckDB taxonomy database. If provided, centrifuger output is annotated with TAX_NAME, KINGDOM, and TAX_RANK.",
            patterns: ["*.duckdb"]
        }
        resolve_strains: {
            description: "When true and taxonomy_db is provided, reclassify 'no rank' nodes below species as 'strain'."
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

    if (defined(taxonomy_db)) {
        call metagenomics.parse_centrifuger_reads as parse_reads {
            input:
                centrifuger_output = run_centrifuger.classification_tsvs[0],
                taxonomy_db        = select_first([taxonomy_db]),
                resolve_strains    = resolve_strains
        }
    }

    output {
        File        centrifuger_reads_classified = select_first([parse_reads.centrifuger_reads_classified, run_centrifuger.classification_tsvs[0]])
    }
}
