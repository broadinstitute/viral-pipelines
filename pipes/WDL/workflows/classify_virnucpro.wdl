version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_virnucpro {
    meta {
        description: "Classifies nucleotide FASTA sequences with VirNucPro and optionally summarizes chunk-level predictions into contig-level calls."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File    in_fasta

        Int     expected_length = 500
        Boolean parallel = false
        Boolean persistent_models = false
        Boolean resume = false
        Boolean v1_fallback = false
        Boolean v1_attention = false
        Int?    batch_size
        Int?    esm_batch_size
        Int?    dnabert_batch_size
        String? gpus

        Boolean run_contig_classification = true
        Float   min_viral_prop = 0.1
        Float   min_nonviral_prop = 0.1
        Int     min_chunks = 5
        String  id_col = "Modified_ID"
        String  id_pattern = "(NODE_[0-9]+)"

        Int     machine_mem_gb = 64
        Int     cpu = 8
        String? accelerator_type
        Int?    accelerator_count
        String? gpu_type
        Int?    gpu_count
        String? vm_size
    }

    String out_basename = basename(basename(basename(in_fasta, ".fasta"), ".fa"), ".fna")

    parameter_meta {
        in_fasta: {
            description: "Input nucleotide sequences in FASTA format.",
            patterns: ["*.fasta", "*.fa", "*.fna"],
            category: "required"
        }
        expected_length: {
            description: "Expected sequence length for the VirNucPro model. Must be 300 or 500.",
            choices: [300, 500],
            category: "common"
        }
        parallel: {
            description: "Enable VirNucPro multi-GPU parallel processing.",
            category: "advanced"
        }
        persistent_models: {
            description: "Keep models resident in GPU memory between stages.",
            category: "advanced"
        }
        resume: {
            description: "Resume from VirNucPro checkpoints.",
            category: "advanced"
        }
        v1_fallback: {
            description: "Use VirNucPro v1.0 multi-worker architecture for ESM-2 instead of v2.0 async DataLoader.",
            category: "advanced"
        }
        v1_attention: {
            description: "Use VirNucPro v1.0-compatible standard attention for ESM-2. Slower, but useful for exact v1 compatibility.",
            category: "advanced"
        }
        batch_size: {
            description: "VirNucPro prediction batch size.",
            category: "advanced"
        }
        esm_batch_size: {
            description: "ESM token batch size.",
            category: "advanced"
        }
        dnabert_batch_size: {
            description: "DNABERT batch size.",
            category: "advanced"
        }
        gpus: {
            description: "Comma-separated GPU IDs to expose to VirNucPro, for example 0,1. Set gpu_count or accelerator_count consistently when using multiple GPUs.",
            category: "advanced"
        }
        run_contig_classification: {
            description: "Summarize chunk-level highest-score predictions into contig-level calls. Disable for non-chunked or non-SPAdes-like FASTA IDs.",
            category: "common"
        }
        min_viral_prop: {
            description: "Minimum confident viral chunk proportion.",
            category: "advanced"
        }
        min_nonviral_prop: {
            description: "Minimum confident non-viral chunk proportion.",
            category: "advanced"
        }
        min_chunks: {
            description: "Minimum chunk count for high/moderate confidence tiers.",
            category: "advanced"
        }
        id_col: {
            description: "Column containing VirNucPro chunk or contig IDs.",
            category: "advanced"
        }
        id_pattern: {
            description: "Regex used to extract contig group IDs from id_col.",
            category: "advanced"
        }
        machine_mem_gb: {
            description: "Memory allocation in GB.",
            category: "runtime"
        }
        cpu: {
            description: "CPU cores to request and pass to VirNucPro.",
            category: "runtime"
        }
        accelerator_type: {
            description: "[GCP/PAPIv2] GPU model to request, for example nvidia-tesla-t4.",
            category: "runtime"
        }
        accelerator_count: {
            description: "[GCP/PAPIv2] Number of GPUs to request.",
            category: "runtime"
        }
        gpu_type: {
            description: "[Terra] GPU model to request, for example nvidia-tesla-t4.",
            category: "runtime"
        }
        gpu_count: {
            description: "[Terra] Number of GPUs to request.",
            category: "runtime"
        }
        vm_size: {
            description: "[TES/Azure] GPU VM size.",
            category: "runtime"
        }
    }

    call metagenomics.virnucpro {
        input:
            in_fasta           = in_fasta,
            expected_length    = expected_length,
            parallel           = parallel,
            persistent_models  = persistent_models,
            resume             = resume,
            v1_fallback        = v1_fallback,
            v1_attention       = v1_attention,
            batch_size         = batch_size,
            esm_batch_size     = esm_batch_size,
            dnabert_batch_size = dnabert_batch_size,
            gpus               = gpus,
            machine_mem_gb     = machine_mem_gb,
            cpu                = cpu,
            accelerator_type   = accelerator_type,
            accelerator_count  = accelerator_count,
            gpu_type           = gpu_type,
            gpu_count          = gpu_count,
            vm_size            = vm_size
    }

    if (run_contig_classification) {
        call metagenomics.virnucpro_contigs {
            input:
                highestscore_tsv = virnucpro.highestscore_tsv,
                out_basename     = out_basename,
                min_viral_prop   = min_viral_prop,
                min_nonviral_prop = min_nonviral_prop,
                min_chunks       = min_chunks,
                id_col           = id_col,
                id_pattern       = id_pattern
        }
    }

    output {
        File   virnucpro_predictions_tsv             = virnucpro.predictions_tsv
        File   virnucpro_highestscore_tsv            = virnucpro.highestscore_tsv
        File?  virnucpro_contig_classifications_tsv  = virnucpro_contigs.contig_classifications_tsv
        Int    virnucpro_max_ram_gb                  = virnucpro.max_ram_gb
        String virnucpro_version                     = virnucpro.virnucpro_version
        String? viralngs_version                     = virnucpro_contigs.viralngs_version
    }
}
