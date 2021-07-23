version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_intrahost.wdl" as intrahost
import "../tasks/tasks_utils.wdl" as utils

workflow mafft_and_snp_annotated {
    meta {
        description: "Align assemblies with mafft and find SNPs with snp-sites."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]     assembly_fastas
        File            ref_fasta
        Int             min_unambig_genome
        Boolean         run_iqtree=false
    }

    parameter_meta {
        assembly_fastas: {
          description: "Set of assembled genomes to align and build trees. These must represent a single chromosome/segment of a genome only. Fastas may be one-sequence-per-individual or a concatenated multi-fasta (unaligned) or a mixture of the two. They may be compressed (gz, bz2, zst, lz4), uncompressed, or a mixture.",
          patterns: ["*.fasta", "*.fa", "*.fasta.gz", "*.fasta.zst"]
        }
        ref_fasta: {
          description: "A reference assembly (not included in assembly_fastas) to align assembly_fastas against. Typically from NCBI RefSeq or similar. Uncompressed.",
          patterns: ["*.fasta", "*.fa"]
        }
        min_unambig_genome: {
          description: "Minimum number of called bases in genome to pass prefilter."
        }
    }

    call utils.zcat {
        input:
            infiles     = assembly_fastas,
            output_name = "all_samples_combined_assembly.fasta.gz"
    }
    call nextstrain.filter_sequences_by_length {
        input:
            sequences_fasta = zcat.combined,
            min_non_N       = min_unambig_genome
    }
    call nextstrain.mafft_one_chr as mafft {
        input:
            sequences = filter_sequences_by_length.filtered_fasta,
            ref_fasta = ref_fasta,
            basename  = "all_samples_aligned.fasta"
    }
    call nextstrain.snp_sites {
        input:
            msa_fasta = mafft.aligned_sequences
    }

    call intrahost.annotate_vcf_snpeff as annotate_vcf {
        input:
            ref_fasta = ref_fasta,
            in_vcf    = snp_sites.snps_vcf
    }

    if(run_iqtree) {
        call nextstrain.draft_augur_tree {
            input:
                msa_or_vcf = mafft.aligned_sequences
        }
    }

    output {
        File  combined_assemblies         = filter_sequences_by_length.filtered_fasta
        File  multiple_alignment          = mafft.aligned_sequences
        File  unmasked_snps               = snp_sites.snps_vcf
        File  unmasked_snps_annotated     = annotate_vcf.annot_vcf_gz
        File  unmasked_snps_annotated_tbi = annotate_vcf.annot_vcf_gz_tbi
        File? ml_tree                     = draft_augur_tree.aligned_tree
    }
}
