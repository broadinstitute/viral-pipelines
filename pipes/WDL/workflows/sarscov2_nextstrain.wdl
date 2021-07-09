version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_utils.wdl" as utils

import "sarscov2_nextstrain_aligned_input.wdl"

workflow sarscov2_nextstrain {
    meta {
        description: "Align assemblies, build trees, and convert to json representation suitable for Nextstrain visualization. See https://nextstrain.org/docs/getting-started/ and https://nextstrain-augur.readthedocs.io/en/stable/"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+    assembly_fastas
        Array[File]+    sample_metadata_tsvs
        File?           ref_fasta
        Int             min_unambig_genome = 27000
    }

    parameter_meta {
        assembly_fastas: {
          description: "Set of assembled genomes to align and build trees. These must represent a single chromosome/segment of a genome only. Fastas may be one-sequence-per-individual or a concatenated multi-fasta (unaligned) or a mixture of the two. They may be compressed (gz, bz2, zst, lz4), uncompressed, or a mixture.",
          patterns: ["*.fasta", "*.fa", "*.fasta.gz", "*.fasta.zst"]
        }
        sample_metadata_tsvs: {
            description: "Tab-separated metadata file that contain binning variables and values. Must contain all samples: output will be filtered to the IDs present in this file.",
            patterns: ["*.txt", "*.tsv"]
        }
        ref_fasta: {
          description: "A reference assembly (not included in assembly_fastas) to align assembly_fastas against. Typically from NCBI RefSeq or similar.",
          patterns: ["*.fasta", "*.fa"]
        }
        min_unambig_genome: {
          description: "Minimum number of called bases in genome to pass prefilter."
        }
    }

    call nextstrain.nextstrain_ncov_defaults

    #### mafft_and_snp

    call utils.gzcat {
        input:
            infiles     = assembly_fastas,
            output_name = "all_samples_combined_assembly.fasta"
    }

    call nextstrain.filter_sequences_by_length {
        input:
            sequences_fasta = gzcat.combined,
            min_non_N       = min_unambig_genome
    }

    call nextstrain.mafft_one_chr_chunked as mafft {
        input:
            sequences = filter_sequences_by_length.filtered_fasta,
            ref_fasta = select_first([ref_fasta, nextstrain_ncov_defaults.reference_fasta]),
            basename  = "all_samples_aligned.fasta"
    }

    call sarscov2_nextstrain_aligned_input.sarscov2_nextstrain_aligned_input {
        input:
            aligned_sequences_fasta = [mafft.aligned_sequences],
            sample_metadata_tsvs    = sample_metadata_tsvs,
            ref_fasta = ref_fasta
    }

    output {
      File             combined_assemblies  = gzcat.combined
      File             multiple_alignment   = mafft.aligned_sequences
      File             unmasked_snps        = sarscov2_nextstrain_aligned_input.unmasked_snps
      
      File             metadata_merged      = sarscov2_nextstrain_aligned_input.metadata_merged
      File             keep_list            = sarscov2_nextstrain_aligned_input.keep_list
      File             subsampled_sequences = sarscov2_nextstrain_aligned_input.subsampled_sequences
      Int              sequences_kept       = sarscov2_nextstrain_aligned_input.sequences_kept
      Map[String, Int] counts_by_group      = sarscov2_nextstrain_aligned_input.counts_by_group
      
      File             ml_tree              = sarscov2_nextstrain_aligned_input.ml_tree
      File             time_tree            = sarscov2_nextstrain_aligned_input.time_tree
      Array[File]      node_data_jsons      = sarscov2_nextstrain_aligned_input.node_data_jsons
      File             tip_frequencies_json = sarscov2_nextstrain_aligned_input.tip_frequencies_json
      File             root_sequence_json   = sarscov2_nextstrain_aligned_input.root_sequence_json
      File             auspice_input_json   = sarscov2_nextstrain_aligned_input.auspice_input_json
    }
}