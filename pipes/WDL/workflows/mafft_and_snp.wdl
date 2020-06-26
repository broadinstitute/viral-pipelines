version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow mafft_and_snp {
    meta {
        description: "Align assemblies with mafft and find SNPs with snp-sites."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]     assembly_fastas
        File            ref_fasta
        Boolean         run_iqtree=false
    }

    parameter_meta {
        assembly_fastas: {
          description: "Set of assembled genomes to align and build trees. These must represent a single chromosome/segment of a genome only. Fastas may be one-sequence-per-individual or a concatenated multi-fasta (unaligned) or a mixture of the two.",
          patterns: ["*.fasta", "*.fa"]
        }
        ref_fasta: {
          description: "A reference assembly (not included in assembly_fastas) to align assembly_fastas against. Typically from NCBI RefSeq or similar.",
          patterns: ["*.fasta", "*.fa"]
        }
    }

    call nextstrain.concatenate {
        input:
            infiles     = assembly_fastas,
            output_name = "all_samples_combined_assembly.fasta"
    }
    call nextstrain.mafft_one_chr as mafft {
        input:
            sequences = concatenate.combined,
            ref_fasta = ref_fasta,
            basename  = "all_samples_aligned.fasta"
    }
    call nextstrain.snp_sites {
        input:
            msa_fasta = mafft.aligned_sequences
    }
    if(run_iqtree) {
        call nextstrain.draft_augur_tree {
            input:
                msa_or_vcf = mafft.aligned_sequences
        }
    }

    output {
        File  combined_assemblies = concatenate.combined
        File  multiple_alignment  = mafft.aligned_sequences
        File  unmasked_snps       = snp_sites.snps_vcf
        File? ml_tree             = draft_augur_tree.aligned_tree
    }
}
