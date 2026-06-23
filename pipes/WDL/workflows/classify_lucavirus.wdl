version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_lucavirus {
    meta {
        description: "Runs LucaVirus protein sequence classification for one input FASTA and one LucaVirus task profile. Input sequences must already be amino-acid/protein sequences; nucleotide contigs should be translated or ORF-called upstream."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File    input_fasta
        String? sample_id
        String  task_profile = "rdrp"
        Int     gpu_id = 0

        Int     helper_machine_mem_gb = 4
        String  viralngs_docker = "quay.io/broadinstitute/viral-ngs:feature-lucavirus-classify"

        Int     machine_mem_gb = 64
        Int     cpu = 8
        String  lucavirus_docker = "ghcr.io/broadinstitute/lucavirus-cuda:v1.0"
        String? accelerator_type
        Int?    accelerator_count
        String? gpu_type
        Int?    gpu_count
        String? predefined_machine_type
        String? vm_size
        Int     boot_disk_size_gb = 100
        Int     disk_size_gb = 50
        Int     preemptible_attempts = 0
    }

    parameter_meta {
        input_fasta: {
            description: "Amino-acid/protein FASTA containing sequences to score with LucaVirus. Nucleotide contigs must be translated or ORF-called upstream.",
            patterns: ["*.fasta", "*.fa", "*.faa"],
            category: "required"
        }
        sample_id: {
            description: "Optional sample identifier used as the output filename prefix. Defaults to the input FASTA basename.",
            category: "common"
        }
        task_profile: {
            description: "Named LucaVirus task profile to run.",
            choices: ["rdrp", "viral_capsid", "virus_ec4"],
            category: "common"
        }
        gpu_id: {
            description: "GPU index to pass to lucavirus-cuda.",
            category: "advanced"
        }
        helper_machine_mem_gb: {
            description: "Memory in GB for viral-ngs helper tasks.",
            category: "runtime"
        }
        viralngs_docker: {
            description: "viral-ngs classify image containing lucavirus_prepare, lucavirus_empty_predictions, and lucavirus_normalize.",
            category: "runtime"
        }
        machine_mem_gb: {
            description: "Memory in GB for the LucaVirus CUDA task.",
            category: "runtime"
        }
        cpu: {
            description: "CPU cores to request for the LucaVirus CUDA task.",
            category: "runtime"
        }
        lucavirus_docker: {
            description: "Standalone CUDA-enabled LucaVirus Docker image with bundled assets.",
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
        predefined_machine_type: {
            description: "[GCP/Terra] Optional exact GCP machine type to request. Required for A100 40GB GPUs; use a2-highgpu-1g for one nvidia-tesla-a100.",
            category: "runtime"
        }
        vm_size: {
            description: "[TES/Azure] GPU VM size.",
            category: "runtime"
        }
        boot_disk_size_gb: {
            description: "Boot disk size in GB for the LucaVirus CUDA task.",
            category: "runtime"
        }
        disk_size_gb: {
            description: "Local working disk size in GB for the LucaVirus CUDA task.",
            category: "runtime"
        }
        preemptible_attempts: {
            description: "Number of preemptible attempts for the LucaVirus CUDA task. Default 0 because LucaVirus has no checkpoint/resume behavior.",
            category: "runtime"
        }
    }

    call metagenomics.lucavirus_prepare {
        input:
            input_fasta    = input_fasta,
            sample_id      = sample_id,
            machine_mem_gb = helper_machine_mem_gb,
            docker         = viralngs_docker
    }

    if (lucavirus_prepare.has_lucavirus_input) {
        call metagenomics.lucavirus {
            input:
                lucavirus_input_csv  = lucavirus_prepare.lucavirus_input_csv,
                output_basename      = lucavirus_prepare.output_basename,
                task_profile         = task_profile,
                gpu_id               = gpu_id,
                machine_mem_gb       = machine_mem_gb,
                cpu                  = cpu,
                docker               = lucavirus_docker,
                accelerator_type     = accelerator_type,
                accelerator_count    = accelerator_count,
                gpu_type             = gpu_type,
                gpu_count            = gpu_count,
                predefined_machine_type = predefined_machine_type,
                vm_size              = vm_size,
                boot_disk_size_gb    = boot_disk_size_gb,
                disk_size_gb         = disk_size_gb,
                preemptible_attempts = preemptible_attempts
        }
    }

    call metagenomics.lucavirus_normalize {
        input:
            raw_predictions_tsv = select_first([
                lucavirus.raw_predictions_tsv,
                lucavirus_prepare.empty_raw_predictions_tsv
            ]),
            output_basename     = lucavirus_prepare.output_basename,
            task_profile        = task_profile,
            machine_mem_gb      = helper_machine_mem_gb,
            docker              = viralngs_docker
    }

    output {
        File    lucavirus_predictions_tsv   = lucavirus_normalize.predictions_tsv
        File    lucavirus_prepare_stats_tsv = lucavirus_prepare.prepare_stats_tsv
        Int     lucavirus_input_sequences   = lucavirus_prepare.n_sequences
        String  viralngs_version            = lucavirus_normalize.viralngs_version
        String? lucavirus_version           = lucavirus.lucavirus_version
        Int?    lucavirus_max_ram_gb        = lucavirus.max_ram_gb
    }
}
