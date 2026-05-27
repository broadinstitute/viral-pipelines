version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_kallisto {
    meta {
        description: "Kallisto/kb pseudoalignment classification of a single BAM or single FASTQ file. Emits long-form counts and read-level Kallisto summary TSV outputs."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File    reads_bam
        File    kallisto_index
        File    t2g
        String? sample_id
        File?   id_to_taxon_map

        Int     kmer_size = 31
        String  parity = "single"
        String  technology = "bulk"
        Boolean protein = false
        Int     extract_threshold = 1
        String  taxonomy_level = "highest"

        Int     machine_mem_gb = 32
        Int     cpu = 16
        String  docker = "quay.io/broadinstitute/viral-ngs:3.0.13-classify"
    }

    parameter_meta {
        reads_bam: {
            description: "Reads to classify and extract from. May be unaligned BAM or single FASTQ/FASTQ.GZ.",
            patterns: ["*.bam", "*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"],
            category: "required"
        }
        kallisto_index: {
            description: "Pre-built Kallisto index file.",
            patterns: ["*.idx", "*.index"],
            category: "required"
        }
        t2g: {
            description: "Transcript-to-gene mapping file. Two-column TSV with transcript IDs in the first column and collapsed target IDs in the second column.",
            patterns: ["*.tsv", "*.txt"],
            category: "required"
        }
        sample_id: {
            description: "Optional sample identifier to stamp into Kallisto count and summary TSV outputs. Defaults to the input reads basename with common BAM/FASTQ extensions removed.",
            category: "common"
        }
        id_to_taxon_map: {
            description: "Optional CSV/TSV mapping Kallisto hit IDs to taxonomy columns. When provided, taxonomy lineage and selected taxonomy name are added to summary.tsv.",
            patterns: ["*.csv", "*.tsv", "*.csv.gz", "*.tsv.gz"],
            category: "common"
        }
        kmer_size: {
            description: "K-mer size used by the Kallisto index. Must match the value used to build the index. Default 31.",
            category: "advanced"
        }
        parity: {
            description: "Library parity passed to kb count. Common values are single or paired. Default single.",
            category: "advanced"
        }
        technology: {
            description: "Technology preset passed to kb count. Default bulk.",
            category: "advanced"
        }
        protein: {
            description: "Indicates that the Kallisto index was built from protein sequences. Default false.",
            category: "advanced"
        }
        extract_threshold: {
            description: "Minimum count threshold used to select targets from the generated Kallisto count h5ad before read extraction. Default 1.",
            category: "advanced"
        }
        taxonomy_level: {
            description: "Taxonomy level to report from id_to_taxon_map. Valid values are highest or deepest. Default highest.",
            category: "advanced"
        }
        machine_mem_gb: {
            description: "Memory allocation in GB for both Kallisto count and extract.",
            category: "runtime"
        }
        cpu: {
            description: "Number of CPUs to request and pass to Kallisto/kb.",
            category: "runtime"
        }
        docker: {
            description: "viral-ngs classify-flavored image containing the refactored `metagenomics kallisto` and `metagenomics kallisto_extract` commands.",
            category: "runtime"
        }
    }

    String derived_sample_id = sub(sub(sub(sub(sub(basename(reads_bam), "\\.bam$", ""), "\\.fastq\\.gz$", ""), "\\.fq\\.gz$", ""), "\\.fastq$", ""), "\\.fq$", "")
    String resolved_sample_id = select_first([sample_id, derived_sample_id])

    call metagenomics.kallisto {
        input:
            reads_bam      = reads_bam,
            kallisto_index = kallisto_index,
            t2g            = t2g,
            sample_id      = resolved_sample_id,
            kmer_size      = kmer_size,
            parity         = parity,
            technology     = technology,
            protein        = protein,
            machine_mem_gb = machine_mem_gb,
            cpu            = cpu,
            docker         = docker
    }

    call metagenomics.kallisto_read_summary {
        input:
            reads_bam      = reads_bam,
            kallisto_index = kallisto_index,
            t2g            = t2g,
            sample_id      = resolved_sample_id,
            h5ad_file      = kallisto.kallisto_counts_h5ad,
            id_to_taxon_map = id_to_taxon_map,
            threshold      = extract_threshold,
            taxonomy_level = taxonomy_level,
            protein        = protein,
            machine_mem_gb = machine_mem_gb,
            cpu            = cpu,
            docker         = docker
    }

    output {
        File kallisto_counts_tsv  = kallisto.kallisto_counts_tsv
        File kallisto_summary_tsv = kallisto_read_summary.summary_tsv
    }
}
