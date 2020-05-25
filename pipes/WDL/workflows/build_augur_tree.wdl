version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow build_augur_tree {
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
    }

    call nextstrain.concatenate {
        input:
            infiles     = assembly_fastas,
            output_name = "all_samples_combined_assembly.fasta"
    }
    call nextstrain.augur_mafft_align {
        input:
            sequences = concatenate.combined,
            ref_fasta = ref_fasta,
            basename  = virus
    }
    call nextstrain.draft_augur_tree {
        input:
            aligned_fasta  = augur_mafft_align.aligned_sequences,
            basename       = virus
    }
    call nextstrain.refine_augur_tree {
        input:
            raw_tree       = draft_augur_tree.aligned_tree,
            aligned_fasta  = augur_mafft_align.aligned_sequences,
            metadata       = sample_metadata,
            basename       = virus
    }
    if(defined(ancestral_traits_to_infer) && length(select_first([ancestral_traits_to_infer,[]]))>0) {
        call nextstrain.ancestral_traits {
            input:
                tree           = refine_augur_tree.tree_refined,
                metadata       = sample_metadata,
                columns        = select_first([ancestral_traits_to_infer,[]]),
                basename       = virus
        }
    }
    call nextstrain.ancestral_tree {
        input:
            refined_tree   = refine_augur_tree.tree_refined,
            aligned_fasta  = augur_mafft_align.aligned_sequences,
            basename       = virus
    }
    call nextstrain.translate_augur_tree {
        input:
            basename       = virus,
            refined_tree   = refine_augur_tree.tree_refined,
            nt_muts        = ancestral_tree.nt_muts_json,
            genbank_gb     = genbank_gb
    }
    call nextstrain.export_auspice_json {
        input:
            tree            = refine_augur_tree.tree_refined,
            sample_metadata = sample_metadata,
            node_data_jsons = select_all([
                                refine_augur_tree.branch_lengths,
                                ancestral_traits.node_data_json,
                                ancestral_tree.nt_muts_json,
                                translate_augur_tree.aa_muts_json])
    }

    output {
        File  combined_assembly_fasta    = concatenate.combined
        File  augur_aligned_fasta        = augur_mafft_align.aligned_sequences
        File  raw_tree                   = draft_augur_tree.aligned_tree
        File  refined_tree               = refine_augur_tree.tree_refined
        File  branch_lengths             = refine_augur_tree.branch_lengths
        File  json_nt_muts               = ancestral_tree.nt_muts_json
        File  ancestral_sequences_fasta  = ancestral_tree.sequences
        File  json_aa_muts               = translate_augur_tree.aa_muts_json
        File? json_ancestral_traits      = ancestral_traits.node_data_json
        File  auspice_input_json         = export_auspice_json.virus_json
    }
}
