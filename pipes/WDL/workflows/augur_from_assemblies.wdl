version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow augur_from_assemblies {
    meta {
        description: "Align assemblies, build trees, and convert to json representation suitable for Nextstrain visualization. See https://nextstrain.org/docs/getting-started/ and https://nextstrain-augur.readthedocs.io/en/stable/"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]     assembly_fastas
        File            sample_metadata
        String          virus
        File            ref_fasta
        File            genbank_gb
        File?           clades_tsv
        Array[String]?  ancestral_traits_to_infer
    }

    parameter_meta {
        assembly_fastas: {
          description: "Set of assembled genomes to align and build trees. These must represent a single chromosome/segment of a genome only. Fastas may be one-sequence-per-individual or a concatenated multi-fasta (unaligned) or a mixture of the two. Fasta header records need to be pipe-delimited (|) for each metadata value.",
          patterns: ["*.fasta", "*.fa"]
        }
        sample_metadata: {
          description: "Metadata in tab-separated text format. See https://nextstrain-augur.readthedocs.io/en/stable/faq/metadata.html for details.",
          patterns: ["*.txt", "*.tsv"]
        }
        virus: {
          description: "A filename-friendly string that is used as a base for output file names."
        }
        ref_fasta: {
          description: "A reference assembly (not included in assembly_fastas) to align assembly_fastas against. Typically from NCBI RefSeq or similar.",
          patterns: ["*.fasta", "*.fa"]
        }
        genbank_gb: {
          description: "A 'genbank' formatted gene annotation file that is used to calculate coding consequences of observed mutations. Must correspond to the same coordinate space as ref_fasta. Typically downloaded from the same NCBI accession number as ref_fasta.",
          patterns: ["*.gb", "*.gbf"]
        }
        ancestral_traits_to_infer: {
          description: "A list of metadata traits to use for ancestral node inference (see https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/traits.html). Multiple traits may be specified; must correspond exactly to column headers in metadata file. Omitting these values will skip ancestral trait inference, and ancestral nodes will not have estimated values for metadata."
        }
        clades_tsv: {
          description: "A TSV file containing clade mutation positions in four columns: [clade  gene    site    alt]; see: https://nextstrain.org/docs/tutorials/defining-clades",
          patterns: ["*.tsv", "*.txt"]
        }
    }

    call nextstrain.concatenate {
        input:
            infiles     = assembly_fastas,
            output_name = "all_samples_combined_assembly.fasta"
    }
    call nextstrain.filter_subsample_sequences {
            input:
                sequences_fasta     = concatenate.combined,
                sample_metadata_tsv = sample_metadata
    }
    call nextstrain.mafft_one_chr as mafft {
        input:
            sequences = filter_subsample_sequences.filtered_fasta,
            ref_fasta = ref_fasta,
            basename  = virus
    }
    call nextstrain.snp_sites {
        input:
            msa_fasta = mafft.aligned_sequences
    }
    call nextstrain.augur_mask_sites {
        input:
            sequences = mafft.aligned_sequences
    }
    call nextstrain.draft_augur_tree {
        input:
            msa_or_vcf = augur_mask_sites.masked_sequences
    }
    call nextstrain.refine_augur_tree {
        input:
            raw_tree    = draft_augur_tree.aligned_tree,
            msa_or_vcf  = augur_mask_sites.masked_sequences,
            metadata    = sample_metadata
    }
    if(defined(ancestral_traits_to_infer) && length(select_first([ancestral_traits_to_infer,[]]))>0) {
        call nextstrain.ancestral_traits {
            input:
                tree           = refine_augur_tree.tree_refined,
                metadata       = sample_metadata,
                columns        = select_first([ancestral_traits_to_infer,[]])
        }
    }
    call nextstrain.ancestral_tree {
        input:
            tree        = refine_augur_tree.tree_refined,
            msa_or_vcf  = augur_mask_sites.masked_sequences
    }
    call nextstrain.translate_augur_tree {
        input:
            tree        = refine_augur_tree.tree_refined,
            nt_muts     = ancestral_tree.nt_muts_json,
            genbank_gb  = genbank_gb
    }
    if(defined(clades_tsv)) {
        call nextstrain.assign_clades_to_nodes {
            input:
                tree_nwk     = refine_augur_tree.tree_refined,
                nt_muts_json = ancestral_tree.nt_muts_json,
                aa_muts_json = translate_augur_tree.aa_muts_json,
                ref_fasta    = ref_fasta,
                clades_tsv   = select_first([clades_tsv])
        }
    }
    call nextstrain.export_auspice_json {
        input:
            tree            = refine_augur_tree.tree_refined,
            sample_metadata = sample_metadata,
            node_data_jsons = select_all([
                                refine_augur_tree.branch_lengths,
                                ancestral_traits.node_data_json,
                                ancestral_tree.nt_muts_json,
                                translate_augur_tree.aa_muts_json,
                                assign_clades_to_nodes.node_clade_data_json])
    }

    output {
        File  combined_assemblies = concatenate.combined
        File  multiple_alignment  = mafft.aligned_sequences
        File  unmasked_snps       = snp_sites.snps_vcf
        File  masked_alignment    = augur_mask_sites.masked_sequences
        File  ml_tree             = draft_augur_tree.aligned_tree
        File  time_tree           = refine_augur_tree.tree_refined
        File  auspice_input_json  = export_auspice_json.virus_json
    }
}
