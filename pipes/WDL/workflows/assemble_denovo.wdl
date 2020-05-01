version 1.0

import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_assembly.wdl" as assembly

workflow assemble_denovo {
  
  input {
    File         reads_unmapped_bam

    Array[File]  deplete_bmtaggerDbs
    Array[File]  deplete_blastDbs
    Array[File]  deplete_bwaDbs

    File?        filter_to_taxon_db

    Boolean      call_isnvs=false
  }

  if(length(deplete_bmtaggerDbs) + length(deplete_blastDbs) + length(deplete_bwaDbs) > 0) {
    call taxon_filter.deplete_taxa {
      input:
        raw_reads_unmapped_bam = reads_unmapped_bam,
        bmtaggerDbs = deplete_bmtaggerDbs,
        blastDbs = deplete_blastDbs,
        bwaDbs = deplete_bwaDbs
    }
  }

  if(defined(filter_to_taxon_db)) {
    call taxon_filter.filter_to_taxon {
      input:
        reads_unmapped_bam = select_first([deplete_taxa.cleaned_bam, reads_unmapped_bam]),
        lastal_db_fasta = select_first([filter_to_taxon_db])
    }
  }

  call read_utils.rmdup_ubam {
    input:
      reads_unmapped_bam = select_first([filter_to_taxon.taxfilt_bam, deplete_taxa.cleaned_bam, reads_unmapped_bam])
  }

  call assembly.assemble {
    input:
      reads_unmapped_bam = rmdup_ubam.dedup_bam
  }

  call assembly.scaffold {
    input:
      contigs_fasta = assemble.contigs_fasta,
      reads_bam = select_first([filter_to_taxon.taxfilt_bam, deplete_taxa.cleaned_bam, reads_unmapped_bam])
  }

  call assembly.refine_2x_and_plot {
    input:
      assembly_fasta = scaffold.scaffold_fasta,
      reads_unmapped_bam = select_first([deplete_taxa.cleaned_bam, reads_unmapped_bam])
  }

  if(call_isnvs) {
    call intrahost.isnvs_per_sample {
        input:
            assembly_fasta = refine_2x_and_plot.final_assembly_fasta,
            mapped_bam = refine_2x_and_plot.aligned_bam
    }
  }

  output {
    File  final_assembly_fasta        = refine_2x_and_plot.final_assembly_fasta
    File  aligned_only_reads_bam      = refine_2x_and_plot.aligned_only_reads_bam
    File  coverage_plot               = refine_2x_and_plot.coverage_plot
    Int   assembly_length             = refine_2x_and_plot.assembly_length
    Int   assembly_length_unambiguous = refine_2x_and_plot.assembly_length_unambiguous
    Int   reads_aligned               = refine_2x_and_plot.reads_aligned
    Float mean_coverage               = refine_2x_and_plot.mean_coverage

    File?  cleaned_bam               = select_first([deplete_taxa.cleaned_bam, reads_unmapped_bam])
    File?  cleaned_fastqc            = deplete_taxa.cleaned_fastqc
    Int?   depletion_read_count_pre  = deplete_taxa.depletion_read_count_pre
    Int?   depletion_read_count_post = deplete_taxa.depletion_read_count_post

    File?  taxfilt_bam               = filter_to_taxon.taxfilt_bam
    File?  taxfilt_fastqc            = filter_to_taxon.taxfilt_fastqc
    Int?   filter_read_count_post    = filter_to_taxon.filter_read_count_post

    File   dedup_bam                 = rmdup_ubam.dedup_bam
    File   dedup_fastqc              = rmdup_ubam.dedup_fastqc
    Int    dedup_read_count_post     = rmdup_ubam.dedup_read_count_post

    File   contigs_fasta             = assemble.contigs_fasta
    File   subsampBam                = assemble.subsampBam
    Int    subsample_read_count      = assemble.subsample_read_count

    File   scaffold_fasta                        = scaffold.scaffold_fasta
    File   intermediate_scaffold_fasta           = scaffold.intermediate_scaffold_fasta
    File   intermediate_gapfill_fasta            = scaffold.intermediate_gapfill_fasta
    Int    assembly_preimpute_length             = scaffold.assembly_preimpute_length
    Int    assembly_preimpute_length_unambiguous = scaffold.assembly_preimpute_length_unambiguous
    String scaffolding_chosen_ref_name           = scaffold.scaffolding_chosen_ref_name
    File   scaffolding_stats                     = scaffold.scaffolding_stats
    File   scaffolding_alt_contigs               = scaffold.scaffolding_alt_contigs

    File?  isnvsFile                   = isnvs_per_sample.isnvsFile

    File aligned_bam                   = refine_2x_and_plot.aligned_bam
    File aligned_only_reads_bam_idx    = refine_2x_and_plot.aligned_only_reads_bam_idx
    File aligned_only_reads_fastqc     = refine_2x_and_plot.aligned_only_reads_fastqc
    File coverage_tsv                  = refine_2x_and_plot.coverage_tsv
    Int  read_pairs_aligned            = refine_2x_and_plot.read_pairs_aligned
    Int  bases_aligned                 = refine_2x_and_plot.bases_aligned

    String? deplete_viral_classify_version  = deplete_taxa.viralngs_version
    String? taxfilt_viral_classify_version  = filter_to_taxon.viralngs_version
    String  assemble_viral_assemble_version = assemble.viralngs_version
    String  scaffold_viral_assemble_version = scaffold.viralngs_version
    String  refine_viral_assemble_version   = refine_2x_and_plot.viralngs_version
    String? isnvs_viral_phylo_version       = isnvs_per_sample.viralngs_version
  }

}
