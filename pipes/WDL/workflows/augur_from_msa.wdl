version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow augur_from_msa {
    meta {
        description: "Build trees, and convert to json representation suitable for Nextstrain visualization. See https://nextstrain.org/docs/getting-started/ and https://nextstrain-augur.readthedocs.io/en/stable/"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        File            msa_or_vcf
        File            sample_metadata
        String          virus
        File            ref_fasta
        File            genbank_gb
        File?           clades_tsv
        Array[String]?  ancestral_traits_to_infer
    }

    parameter_meta {
        msa_or_vcf: {
          description: "Multiple sequence alignment (aligned fasta) or variants (vcf format).",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
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

    call nextstrain.augur_mask_sites {
        input:
            sequences = msa_or_vcf
    }
    call nextstrain.draft_augur_tree {
        input:
            msa_or_vcf = augur_mask_sites.masked_sequences,
            basename   = virus
    }
    call nextstrain.refine_augur_tree {
        input:
            raw_tree    = draft_augur_tree.aligned_tree,
            msa_or_vcf  = augur_mask_sites.masked_sequences,
            metadata    = sample_metadata,
            basename    = virus
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
            refined_tree  = refine_augur_tree.tree_refined,
            msa_or_vcf    = augur_mask_sites.masked_sequences,
            basename      = virus
    }
    call nextstrain.translate_augur_tree {
        input:
            basename       = virus,
            refined_tree   = refine_augur_tree.tree_refined,
            nt_muts        = ancestral_tree.nt_muts_json,
            genbank_gb     = genbank_gb
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
        File  masked_fasta               = augur_mask_sites.masked_sequences
        File  raw_tree                   = draft_augur_tree.aligned_tree
        File  refined_tree               = refine_augur_tree.tree_refined
        File  branch_lengths             = refine_augur_tree.branch_lengths
        File  json_nt_muts               = ancestral_tree.nt_muts_json
        File  ancestral_sequences_fasta  = ancestral_tree.sequences
        File  json_aa_muts               = translate_augur_tree.aa_muts_json
        File? node_clade_data_json       = assign_clades_to_nodes.node_clade_data_json
        File? json_ancestral_traits      = ancestral_traits.node_data_json
        File  auspice_input_json         = export_auspice_json.virus_json
    }
}
