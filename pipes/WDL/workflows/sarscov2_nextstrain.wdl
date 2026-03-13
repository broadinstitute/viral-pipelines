version 1.1
#DX_SKIP_WORKFLOW

import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_utils.wdl" as utils

import "sarscov2_nextstrain_aligned_input.wdl"

workflow sarscov2_nextstrain {
    meta {
        description: "Align assemblies, build trees, and convert to json representation suitable for Nextstrain visualization. See https://nextstrain.org/docs/getting-started/ and https://nextstrain-augur.readthedocs.io/en/stable/"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]+    assembly_fastas=["gs://nextstrain-data/files/ncov/open/sequences.fasta.zst"]
        Array[File]+    sample_metadata_tsvs=["gs://nextstrain-data/files/ncov/open/metadata.tsv.gz"]
        File?           ref_fasta
        Int             min_unambig_genome = 27000
        String          tree_root_seq_id = "Wuhan-Hu-1/2019"


        String          build_name
        File            builds_yaml

        Array[String]?  ancestral_traits_to_infer

        File?           auspice_config
        File?           clades_tsv
        File?           lat_longs_tsv
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

    call utils.zcat {
        input:
            infiles     = assembly_fastas,
            output_name = "all_samples_combined_assembly.fasta.zst"
    }

    call nextstrain.nextstrain_deduplicate_sequences as dedup_seqs {
        input:
            sequences_fasta = zcat.combined
    }

    call utils.filter_sequences_by_length {
        input:
            sequences_fasta = dedup_seqs.sequences_deduplicated_fasta,
            min_non_N       = min_unambig_genome
    }

    call nextstrain.mafft_one_chr_chunked as mafft {
        input:
            sequences = filter_sequences_by_length.filtered_fasta,
            ref_fasta = select_first([ref_fasta, nextstrain_ncov_defaults.reference_fasta]),
            basename  = "all_samples_aligned.fasta"
    }

    #### merge metadata, compute derived cols
    if(length(sample_metadata_tsvs)>1) {
        call utils.tsv_join {
            input:
                input_tsvs     = sample_metadata_tsvs,
                id_col         = 'strain',
                out_basename   = "metadata-merged",
                machine_mem_gb = 30
        }
    }
    call nextstrain.derived_cols {
        input:
            metadata_tsv = select_first(flatten([[tsv_join.out_tsv], sample_metadata_tsvs]))
    }


    #### subsample sequences with nextstrain yaml file
    call nextstrain.nextstrain_build_subsample as subsample {
        input:
            alignment_msa_fasta = mafft.aligned_sequences,
            sample_metadata_tsv = derived_cols.derived_metadata,
            build_name          = build_name,
            builds_yaml         = builds_yaml
    }
    call utils.fasta_to_ids {
        input:
            sequences_fasta = subsample.subsampled_msa
    }
    call nextstrain.snp_sites {
        input:
            msa_fasta = subsample.subsampled_msa
    }


    #### augur_from_msa

    call nextstrain.augur_mask_sites {
        input:
            sequences = subsample.subsampled_msa
    }

    call nextstrain.draft_augur_tree {
        input:
            msa_or_vcf = augur_mask_sites.masked_sequences
    }

    call nextstrain.refine_augur_tree {
        input:
            raw_tree   = draft_augur_tree.aligned_tree,
            msa_or_vcf = augur_mask_sites.masked_sequences,
            metadata   = derived_cols.derived_metadata,
            root       = tree_root_seq_id
    }

    if(defined(ancestral_traits_to_infer) && length(select_first([ancestral_traits_to_infer,[]]))>0) {
        call nextstrain.ancestral_traits {
            input:
                tree     = refine_augur_tree.tree_refined,
                metadata = derived_cols.derived_metadata,
                columns  = select_first([ancestral_traits_to_infer,[]])
        }
    }

    call nextstrain.tip_frequencies {
        input:
            tree                 = refine_augur_tree.tree_refined,
            metadata             = derived_cols.derived_metadata,
            min_date             = 2020.0,
            pivot_interval       = 1,
            pivot_interval_units = "weeks",
            narrow_bandwidth     = 0.05,
            proportion_wide      = 0.0,
            out_basename         = "auspice-~{build_name}"
    }

    call nextstrain.ancestral_tree {
        input:
            tree       = refine_augur_tree.tree_refined,
            msa_or_vcf = augur_mask_sites.masked_sequences
    }

    call nextstrain.translate_augur_tree {
        input:
            tree       = refine_augur_tree.tree_refined,
            nt_muts    = ancestral_tree.nt_muts_json,
            genbank_gb = nextstrain_ncov_defaults.reference_gb
    }

    call nextstrain.assign_clades_to_nodes {
        input:
            tree_nwk     = refine_augur_tree.tree_refined,
            nt_muts_json = ancestral_tree.nt_muts_json,
            aa_muts_json = translate_augur_tree.aa_muts_json,
            ref_fasta    = select_first([ref_fasta, nextstrain_ncov_defaults.reference_fasta]),
            clades_tsv   = select_first([clades_tsv, nextstrain_ncov_defaults.clades_tsv])
    }

    call nextstrain.export_auspice_json {
        input:
            tree            = refine_augur_tree.tree_refined,
            sample_metadata = derived_cols.derived_metadata,
            lat_longs_tsv   = select_first([lat_longs_tsv, nextstrain_ncov_defaults.lat_longs_tsv]),
            node_data_jsons = select_all([
                                refine_augur_tree.branch_lengths,
                                ancestral_traits.node_data_json,
                                ancestral_tree.nt_muts_json,
                                translate_augur_tree.aa_muts_json,
                                assign_clades_to_nodes.node_clade_data_json]),
            auspice_config  = select_first([auspice_config, nextstrain_ncov_defaults.auspice_config]),
            out_basename    = "auspice-~{build_name}"
    }


    output {
      File             combined_assemblies  = zcat.combined
      File             multiple_alignment   = mafft.aligned_sequences
      File             unmasked_snps        = snp_sites.snps_vcf
      
      File             metadata_merged      = derived_cols.derived_metadata
      File             keep_list            = fasta_to_ids.ids_txt
      File             subsampled_sequences = subsample.subsampled_msa
      File             masked_alignment     = augur_mask_sites.masked_sequences
      Int              sequences_kept       = subsample.sequences_out
      Map[String, Int] counts_by_group      = subsample.counts_by_group
      
      File             ml_tree              = draft_augur_tree.aligned_tree
      File             time_tree            = refine_augur_tree.tree_refined
      Array[File]      node_data_jsons      = select_all([
                                                 refine_augur_tree.branch_lengths,
                                                 ancestral_traits.node_data_json,
                                                 ancestral_tree.nt_muts_json,
                                                 translate_augur_tree.aa_muts_json,
                                                 assign_clades_to_nodes.node_clade_data_json])
      File             tip_frequencies_json = tip_frequencies.node_data_json
      File             root_sequence_json   = export_auspice_json.root_sequence_json
      File             auspice_input_json   = export_auspice_json.virus_json
    }
}
