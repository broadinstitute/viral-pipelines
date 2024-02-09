version 1.0

import "../tasks/tasks_assembly.wdl" as assembly
import "assemble_refbased.wdl" as assemble_refbased

workflow scaffold_and_refine {
    meta {
        description: "Scaffold de novo contigs against a set of possible references and subsequently polish with reads."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File reads_unmapped_bam
    }

    call assembly.scaffold {
        input:
            reads_bam = reads_unmapped_bam
    }

    call assemble_refbased.assemble_refbased as refine {
        input:
            reads_unmapped_bams = [reads_unmapped_bam],
            reference_fasta     = scaffold.scaffold_fasta
    }

  output {
    File   final_assembly_fasta                  = refine.assembly_fasta
    File   aligned_only_reads_bam                = refine.align_to_self_merged_aligned_only_bam
    File   coverage_plot                         = refine.align_to_self_merged_coverage_plot
    Int    assembly_length                       = refine.assembly_length
    Int    assembly_length_unambiguous           = refine.assembly_length_unambiguous
    Int    reads_aligned                         = refine.align_to_self_merged_reads_aligned
    Float  mean_coverage                         = refine.align_to_self_merged_mean_coverage
    
    File   scaffold_fasta                        = scaffold.scaffold_fasta
    File   intermediate_scaffold_fasta           = scaffold.intermediate_scaffold_fasta
    File   intermediate_gapfill_fasta            = scaffold.intermediate_gapfill_fasta
    Int    assembly_preimpute_length             = scaffold.assembly_preimpute_length
    Int    assembly_preimpute_length_unambiguous = scaffold.assembly_preimpute_length_unambiguous
    Array[String]  scaffolding_chosen_ref_names  = scaffold.scaffolding_chosen_ref_names
    File   scaffolding_stats                     = scaffold.scaffolding_stats
    File   scaffolding_alt_contigs               = scaffold.scaffolding_alt_contigs

    Int    replicate_concordant_sites            = refine.replicate_concordant_sites
    Int    replicate_discordant_snps             = refine.replicate_discordant_snps
    Int    replicate_discordant_indels           = refine.replicate_discordant_indels
    Int    num_read_groups                       = refine.num_read_groups
    Int    num_libraries                         = refine.num_libraries
    File   replicate_discordant_vcf              = refine.replicate_discordant_vcf

    File   isnvsFile                             = refine.align_to_self_isnvs_vcf
    
    File   aligned_bam                           = refine.align_to_self_merged_aligned_only_bam
    File   aligned_only_reads_fastqc             = refine.align_to_ref_fastqc
    File   coverage_tsv                          = refine.align_to_self_merged_coverage_tsv
    Int    read_pairs_aligned                    = refine.align_to_self_merged_read_pairs_aligned
    Float  bases_aligned                         = refine.align_to_self_merged_bases_aligned
    
    String scaffold_viral_assemble_version       = scaffold.viralngs_version
    String refine_viral_assemble_version         = refine.viral_assemble_version
  }
}
