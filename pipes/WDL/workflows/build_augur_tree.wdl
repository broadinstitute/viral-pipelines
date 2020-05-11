version 1.0

import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow build_augur_tree {
    input {
        Array[File]     assembly_fastas # fasta header records need to be "|" delimited for each metadata value
        File            metadata
        String          virus
        File            ref_fasta       # reference genome (often RefSeq)
        File            genbank_gb      # Genbank file for amino acid annotations (same coord space as ref_fasta, typically RefSeq)
        Array[String]?  ancestral_traits_to_infer
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
            metadata       = metadata,
            basename       = virus
    }
    if(defined(ancestral_traits_to_infer) && length(select_first([ancestral_traits_to_infer,[]]))>0) {
        call nextstrain.ancestral_traits {
            input:
                tree           = refine_augur_tree.tree_refined,
                metadata       = metadata,
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
            refined_tree   = refine_augur_tree.tree_refined,
            metadata       = metadata,
            branch_lengths = refine_augur_tree.branch_lengths,
            traits         = ancestral_traits.node_data_json,
            nt_muts        = ancestral_tree.nt_muts_json,
            aa_muts        = translate_augur_tree.aa_muts_json,
            basename       = virus
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
