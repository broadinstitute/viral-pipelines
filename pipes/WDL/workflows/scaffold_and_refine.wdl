version 1.0

import "../tasks/tasks_assembly.wdl" as assembly

workflow scaffold_and_refine {
    input {
        File reads_unmapped_bam
    }

    call assembly.scaffold {
        input:
            reads_bam = reads_unmapped_bam
    }

    call assembly.refine_2x_and_plot {
        input:
            assembly_fasta = scaffold.scaffold_fasta,
            reads_unmapped_bam = reads_unmapped_bam
    }

  output {
    File  final_assembly_fasta        = refine_2x_and_plot.final_assembly_fasta
    File  aligned_only_reads_bam      = refine_2x_and_plot.aligned_only_reads_bam
    File  coverage_plot               = refine_2x_and_plot.coverage_plot
    Int   assembly_length             = refine_2x_and_plot.assembly_length
    Int   assembly_length_unambiguous = refine_2x_and_plot.assembly_length_unambiguous
    Int   reads_aligned               = refine_2x_and_plot.reads_aligned
    Float mean_coverage               = refine_2x_and_plot.mean_coverage

    File   scaffold_fasta                        = scaffold.scaffold_fasta
    File   intermediate_scaffold_fasta           = scaffold.intermediate_scaffold_fasta
    File   intermediate_gapfill_fasta            = scaffold.intermediate_gapfill_fasta
    Int    assembly_preimpute_length             = scaffold.assembly_preimpute_length
    Int    assembly_preimpute_length_unambiguous = scaffold.assembly_preimpute_length_unambiguous
    String scaffolding_chosen_ref_name           = scaffold.scaffolding_chosen_ref_name
    File   scaffolding_stats                     = scaffold.scaffolding_stats
    File   scaffolding_alt_contigs               = scaffold.scaffolding_alt_contigs

    File aligned_bam                   = refine_2x_and_plot.aligned_bam
    File aligned_only_reads_bam_idx    = refine_2x_and_plot.aligned_only_reads_bam_idx
    File aligned_only_reads_fastqc     = refine_2x_and_plot.aligned_only_reads_fastqc
    File coverage_tsv                  = refine_2x_and_plot.coverage_tsv
    Int  read_pairs_aligned            = refine_2x_and_plot.read_pairs_aligned
    Float bases_aligned                 = refine_2x_and_plot.bases_aligned

    String scaffold_viral_assemble_version = scaffold.viralngs_version
    String refine_viral_assemble_version   = refine_2x_and_plot.viralngs_version
  }
}
