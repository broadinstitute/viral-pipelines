version 1.0

import "mafft_and_snp.wdl"
import "subsample_by_metadata_with_focal.wdl"
import "augur_from_msa.wdl"


workflow augur_from_assemblies {
    meta {
        description: "Align assemblies, build trees, and convert to json representation suitable for Nextstrain visualization. See https://nextstrain.org/docs/getting-started/ and https://nextstrain-augur.readthedocs.io/en/stable/"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        File            ref_fasta
    }

    parameter_meta {
        ref_fasta: {
          description: "A reference assembly (not included in assembly_fastas) to align assembly_fastas against. Typically from NCBI RefSeq or similar.",
          patterns: ["*.fasta", "*.fa"]
        }
    }

    call mafft_and_snp.mafft_and_snp {
        input:
            ref_fasta = ref_fasta
    }

    call subsample_by_metadata_with_focal.subsample_by_metadata_with_focal {
        input:
            sequences_fasta = mafft_and_snp.multiple_alignment
    }

    call augur_from_msa.augur_from_msa {
        input:
            msa_or_vcf = mafft_and_snp.multiple_alignment,
            sample_metadata = [subsample_by_metadata_with_focal.metadata_merged],
            ref_fasta = ref_fasta,
            keep_list = [subsample_by_metadata_with_focal.keep_list]
    }

    output {
        File  combined_assemblies  = mafft_and_snp.combined_assemblies
        File  multiple_alignment   = mafft_and_snp.multiple_alignment
        File  unmasked_snps        = mafft_and_snp.unmasked_snps

        File  metadata_merged      = subsample_by_metadata_with_focal.metadata_merged
        File  keep_list            = subsample_by_metadata_with_focal.keep_list
        File  subsampled_sequences = subsample_by_metadata_with_focal.subsampled_sequences
        Int   focal_kept           = subsample_by_metadata_with_focal.focal_kept
        Int   global_kept          = subsample_by_metadata_with_focal.global_kept
        Int   sequences_kept       = subsample_by_metadata_with_focal.sequences_kept

        File  masked_alignment     = augur_from_msa.masked_alignment
        File  ml_tree              = augur_from_msa.ml_tree
        File  time_tree            = augur_from_msa.time_tree
        File  auspice_input_json   = augur_from_msa.auspice_input_json
    }
}