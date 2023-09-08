version 1.0

import "../tasks/tasks_interhost.wdl" as interhost
import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_utils.wdl" as utils

workflow augur_from_msa_with_subsampler {
    meta {
        description: "Build trees, and convert to json representation suitable for Nextstrain visualization. See https://nextstrain.org/docs/getting-started/ and https://nextstrain-augur.readthedocs.io/en/stable/"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File           aligned_msa_fasta
        Array[File]+   sample_metadata
        File?          ref_fasta
        File?          genbank_gb
        File           auspice_config
        File?          clades_tsv
        Array[String]? ancestral_traits_to_infer
        File?          mask_bed
    }

    parameter_meta {
        aligned_msa_fasta: {
          description: "Multiple sequence alignment (aligned fasta).",
          patterns: ["*.fasta", "*.fa"]
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
        mask_bed: {
          description: "Optional list of sites to mask when building trees.",
          patterns: ["*.bed"]
        }
    }

    # merge tsvs if necessary
    if(length(sample_metadata)>1) {
        call utils.tsv_join {
            input:
                input_tsvs   = sample_metadata,
                id_col       = 'strain',
                out_basename = "metadata-merged"
        }
    }

    # subsample and filter genomic data based on epi case data
    call interhost.subsample_by_cases {
        input:
            metadata = select_first(flatten([[tsv_join.out_tsv], sample_metadata]))
    }
    call nextstrain.filter_sequences_to_list {
        input:
            sequences = aligned_msa_fasta,
            keep_list = [subsample_by_cases.selected_sequences]
    }

    # standard augur pipeline
    if(defined(mask_bed)) {
        call nextstrain.augur_mask_sites {
            input:
                sequences = filter_sequences_to_list.filtered_fasta,
                mask_bed  = mask_bed
        }
    }
    File masked_sequences = select_first([augur_mask_sites.masked_sequences, filter_sequences_to_list.filtered_fasta])
    call nextstrain.draft_augur_tree {
        input:
            msa_or_vcf = masked_sequences
    }
    call nextstrain.refine_augur_tree {
        input:
            raw_tree   = draft_augur_tree.aligned_tree,
            msa_or_vcf = masked_sequences,
            metadata   = subsample_by_cases.selected_metadata
    }
    if(defined(ancestral_traits_to_infer) && length(select_first([ancestral_traits_to_infer,[]]))>0) {
        call nextstrain.ancestral_traits {
            input:
                tree     = refine_augur_tree.tree_refined,
                metadata = subsample_by_cases.selected_metadata,
                columns  = select_first([ancestral_traits_to_infer,[]])
        }
    }
    call nextstrain.tip_frequencies {
        input:
            tree     = refine_augur_tree.tree_refined,
            metadata = subsample_by_cases.selected_metadata
    }
    call nextstrain.ancestral_tree {
        input:
            tree       = refine_augur_tree.tree_refined,
            msa_or_vcf = masked_sequences
    }
    if(defined(genbank_gb)) {
        call nextstrain.translate_augur_tree {
            input:
                tree       = refine_augur_tree.tree_refined,
                nt_muts    = ancestral_tree.nt_muts_json,
                genbank_gb = select_first([genbank_gb])
        }
    }
    if(defined(clades_tsv) && defined(ref_fasta)) {
        call nextstrain.assign_clades_to_nodes {
            input:
                tree_nwk     = refine_augur_tree.tree_refined,
                nt_muts_json = ancestral_tree.nt_muts_json,
                aa_muts_json = translate_augur_tree.aa_muts_json,
                ref_fasta    = select_first([ref_fasta]),
                clades_tsv   = select_first([clades_tsv])
        }
    }
    call nextstrain.export_auspice_json {
        input:
            tree            = refine_augur_tree.tree_refined,
            sample_metadata = subsample_by_cases.selected_metadata,
            node_data_jsons = select_all([
                                refine_augur_tree.branch_lengths,
                                ancestral_traits.node_data_json,
                                ancestral_tree.nt_muts_json,
                                translate_augur_tree.aa_muts_json,
                                assign_clades_to_nodes.node_clade_data_json]),
            auspice_config  = auspice_config
    }

    output {
        File        masked_subsampled_msa = masked_sequences
        
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