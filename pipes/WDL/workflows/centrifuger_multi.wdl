version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow centrifuger_multi {
    meta {
        description: "Batch Centrifuger taxonomic classification of multiple BAM files. All samples are classified on a single node with the database loaded once, amortizing the 200+ GB index load across all samples."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]+  reads_bams
        File          centrifuger_db_tgz
        String        db_name

        File?         taxonomy_db
        Boolean       resolve_strains = false

        Int           machine_mem_gb = 256
        Int           cpu            = 16
        String        docker = "ghcr.io/broadinstitute/docker-centrifuger:1.0.0"
    }

    parameter_meta {
        reads_bams: {
            description: "Array of BAM files, one per sample. Each BAM is converted to FASTQ and classified sequentially on a single node.",
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
            description: "Pre-built DuckDB taxonomy database. If provided, each centrifuger output is annotated with TAX_NAME, KINGDOM, and TAX_RANK.",
            patterns: ["*.duckdb"]
        }
        resolve_strains: {
            description: "When true and taxonomy_db is provided, reclassify 'no rank' nodes below species as 'strain'."
        }
        machine_mem_gb: {
            description: "Memory in GB. Default 256 GB sized for multi-sample batch with NT-scale index."
        }
        cpu: {
            description: "Number of CPUs. Default 16 for multi-sample batch processing."
        }
    }

    call metagenomics.centrifuger as run_centrifuger {
        input:
            reads_bams         = reads_bams,
            centrifuger_db_tgz = centrifuger_db_tgz,
            db_name            = db_name,
            machine_mem_gb     = machine_mem_gb,
            cpu                = cpu,
            docker             = docker
    }

    if (defined(taxonomy_db)) {
        scatter (tsv in run_centrifuger.classification_tsvs) {
            call metagenomics.parse_centrifuger_reads as parse_reads {
                input:
                    centrifuger_output = tsv,
                    taxonomy_db        = select_first([taxonomy_db]),
                    resolve_strains    = resolve_strains
            }
        }
    }

    output {
        Array[File]  centrifuger_reads_classified = select_first([parse_reads.centrifuger_reads_classified, run_centrifuger.classification_tsvs])
    }
}
