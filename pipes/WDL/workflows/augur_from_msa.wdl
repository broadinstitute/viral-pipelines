version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_utils.wdl" as utils

workflow augur_from_msa {
    meta {
        description: "Build trees, and convert to json representation suitable for Nextstrain visualization. See https://nextstrain.org/docs/getting-started/ and https://nextstrain-augur.readthedocs.io/en/stable/"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        File           msa_or_vcf
        Array[File]+   sample_metadata
        File           ref_fasta
        File           genbank_gb
        File           auspice_config
        File?          clades_tsv
        Array[String]? ancestral_traits_to_infer
        Array[File]?   keep_list
        File?          mask_bed
    }

    parameter_meta {
        msa_or_vcf: {
          description: "Multiple sequence alignment (aligned fasta) or variants (vcf format).",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
        sample_metadata: {
          description: "Metadata in tab-separated text format. See https://nextstrain-augur.readthedocs.io/en/stable/faq/metadata.html for details. At least one tab file must be provided--if multiple are provided, they will be joined via a full left outer join using the 'strain' column as the join ID.",
          patterns: ["*.txt", "*.tsv"]
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
        auspice_config: {
          description: "A file specifying options to customize the auspice export; see: https://nextstrain.github.io/auspice/customise-client/introduction",
          patterns: ["*.json", "*.txt"]
        }
        clades_tsv: {
          description: "A TSV file containing clade mutation positions in four columns: [clade  gene    site    alt]; see: https://nextstrain.org/docs/tutorials/defining-clades",
          patterns: ["*.tsv", "*.txt"]
        }
        keep_list: {
          description: "Optional lists of strain ids to filter inputs down to.",
          patterns: ["*.txt", "*.tsv"]
        }
        mask_bed: {
          description: "Optional list of sites to mask when building trees.",
          patterns: ["*.bed"]
        }
    }

    call nextstrain.filter_sequences_to_list {
        input:
            sequences = msa_or_vcf,
            keep_list = keep_list
    }
    call nextstrain.augur_mask_sites {
        input:
            sequences = filter_sequences_to_list.filtered_fasta,
            mask_bed  = mask_bed
    }
    call nextstrain.draft_augur_tree {
        input:
            msa_or_vcf = augur_mask_sites.masked_sequences
    }
    if(length(sample_metadata)>1) {
        call utils.tsv_join {
            input:
                input_tsvs   = sample_metadata,
                id_col       = 'strain',
                out_basename = "metadata-merged"
        }
    }
    call nextstrain.refine_augur_tree {
        input:
            raw_tree   = draft_augur_tree.aligned_tree,
            msa_or_vcf = augur_mask_sites.masked_sequences,
            metadata   = select_first(flatten([[tsv_join.out_tsv], sample_metadata]))
    }
    if(defined(ancestral_traits_to_infer) && length(select_first([ancestral_traits_to_infer,[]]))>0) {
        call nextstrain.ancestral_traits {
            input:
                tree     = refine_augur_tree.tree_refined,
                metadata = select_first(flatten([[tsv_join.out_tsv], sample_metadata])),
                columns  = select_first([ancestral_traits_to_infer,[]])
        }
    }
    call nextstrain.tip_frequencies {
        input:
            tree     = refine_augur_tree.tree_refined,
            metadata = select_first(flatten([[tsv_join.out_tsv], sample_metadata]))
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
            genbank_gb = genbank_gb
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
            sample_metadata = select_first(flatten([[tsv_join.out_tsv], sample_metadata])),
            node_data_jsons = select_all([
                                refine_augur_tree.branch_lengths,
                                ancestral_traits.node_data_json,
                                ancestral_tree.nt_muts_json,
                                translate_augur_tree.aa_muts_json,
                                assign_clades_to_nodes.node_clade_data_json]),
            auspice_config  = auspice_config
    }

    output {
        File        masked_alignment      = augur_mask_sites.masked_sequences
        
        File        ml_tree               = draft_augur_tree.aligned_tree
        File        time_tree             = refine_augur_tree.tree_refined
        
        Array[File] node_data_jsons       = select_all([
                    refine_augur_tree.branch_lengths,
                    ancestral_traits.node_data_json,
                    ancestral_tree.nt_muts_json,
                    translate_augur_tree.aa_muts_json,
                    assign_clades_to_nodes.node_clade_data_json])

        File        auspice_input_json    = export_auspice_json.virus_json
        File        tip_frequencies_json  = tip_frequencies.node_data_json
        File        root_sequence_json    = export_auspice_json.root_sequence_json
    }
}