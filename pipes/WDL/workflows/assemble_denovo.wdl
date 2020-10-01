version 1.0

import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_intrahost.wdl" as intrahost

workflow assemble_denovo {
  
  input {
    Array[File]+ reads_unmapped_bams

    Array[File]+ reference_genome_fasta

    Array[File]  deplete_bmtaggerDbs = []
    Array[File]  deplete_blastDbs = []
    Array[File]  deplete_bwaDbs =[]

    File?        filter_to_taxon_db
    File         trim_clip_db

    File?        novocraft_license

    Boolean      call_isnvs=false

    String       assembler="spades"
    Float?       scaffold_min_length_fraction
    Float?       scaffold_min_unambig
    Int?         scaffold_replace_length=55
    Int?         nucmer_max_gap
    Int?         nucmer_min_match
    Int?         nucmer_min_cluster
    Float?       scaffold_min_pct_contig_aligned
  }

  parameter_meta {
    reads_unmapped_bams: {
       description: "Unaligned reads in BAM format",
       patterns: ["*.bam"]
    }
    deplete_bmtaggerDbs: {
       description: "Optional list of databases to use for bmtagger-based depletion. Sequences in fasta format will be indexed on the fly, pre-bmtagger-indexed databases may be provided as tarballs.",
       patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    deplete_blastDbs: {
      description: "Optional list of databases to use for blastn-based depletion. Sequences in fasta format will be indexed on the fly, pre-blast-indexed databases may be provided as tarballs.",
      patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    deplete_bwaDbs: {
      description: "Optional list of databases to use for bwa mem-based depletion. Sequences in fasta format will be indexed on the fly, pre-bwa-indexed databases may be provided as tarballs.",
      patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    filter_to_taxon_db: {
      description: "Optional database to use to filter read set to those that match by LASTAL. Sequences in fasta format will be indexed on the fly.",
      patterns: ["*.fasta"]
    }
    reference_genome_fasta: {
      description: "After denovo assembly, large contigs are scaffolded against a reference genome to determine orientation and to join contigs together, before further polishing by reads. You must supply at least one reference genome (all segments/chromomes in a single fasta file). If more than one reference is provided, contigs will be scaffolded against all of them and the one with the most complete assembly will be chosen for downstream polishing.",
      patterns: ["*.fasta"]
    }
  }

  String sample_name = basename(basename(reads_unmapped_bam, ".bam"), ".cleaned")

  scatter(reads_unmapped_bam in reads_unmapped_bams){
      if(length(deplete_bmtaggerDbs) + length(deplete_blastDbs) + length(deplete_bwaDbs) > 0) {
        call taxon_filter.deplete_taxa {
          input:
            raw_reads_unmapped_bam = reads_unmapped_bam,
            bmtaggerDbs            = deplete_bmtaggerDbs,
            blastDbs               = deplete_blastDbs,
            bwaDbs                 = deplete_bwaDbs
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
      Map[File,File] preprocessed_bams_map = {"ubam":reads_unmapped_bam,"cleaned_bam":select_first([deplete_taxa.cleaned_bam,[]])}
  }

  call read_utils.merge_and_reheader_bams as merge_depleted_taxfilt_rmdup {
      input:
          in_bams      = rmdup_ubam.dedup_bam,
          sample_name  = sample_name,
          out_basename = "${sample_name}.depleted.taxfilt.rmdup"
  }

  call assembly.assemble {
    input:
      reads_unmapped_bam = merge_depleted_taxfilt_rmdup.out_bam,
      trim_clip_db       = trim_clip_db,
      always_succeed     = true,
      assembler          = assembler,
      sample_name        = sample_name
  }

  call assembly.scaffold {
    input:
      contigs_fasta                   = assemble.contigs_fasta,
      reads_bam                       = select_first([filter_to_taxon.taxfilt_bam, deplete_taxa.cleaned_bam, reads_unmapped_bam]),
      reference_genome_fasta          = reference_genome_fasta,
      min_length_fraction             = scaffold_min_length_fraction,
      min_unambig                     = scaffold_min_unambig,
      replace_length                  = scaffold_replace_length,
      nucmer_max_gap                  = nucmer_max_gap,
      nucmer_min_match                = nucmer_min_match,
      nucmer_min_cluster              = nucmer_min_cluster,
      scaffold_min_pct_contig_aligned = scaffold_min_pct_contig_aligned
  }

  # deplete_taxa is stuck in scatter
  # we go through the trouble of this scatter because we need to
  # go back to raw or taxon-depleted reads, and the 
  # previously merged bam is the result of filtering to 
  # like-taxa. Going back to raw reads allows divergent reads
  # to map to the denovo scaffols
  scatter(preprocessed_bams in preprocessed_bams_map){
      call assembly.refine_2x_and_plot {
        input:
          assembly_fasta     = scaffold.scaffold_fasta,
          reads_unmapped_bam = select_first([preprocessed_bams["cleaned_bam"], preprocessed_bams["ubam"]]),
          novocraft_license  = novocraft_license,
          sample_name        = sample_name
      }
  }

  # after going back to starting reads, merge the ones that map to the scaffold
  # before calling refine yet again
  call read_utils.merge_and_reheader_bams as merge_refined {
      input:
          in_bams      = refine_2x_and_plot.aligned_bam,
          sample_name  = sample_name,
          out_basename = "${sample_name}.refine-mapped"
  }

  call assembly.refine_2x_and_plot as final_refinement {
    input:
      assembly_fasta     = scaffold.scaffold_fasta,
      reads_unmapped_bam = merge_refined.out_bam
      novocraft_license  = novocraft_license,
      sample_name        = sample_name
  }

  call assembly.run_discordance {
        input:
            reads_aligned_bam = final_refinement.aligned_bam,
            reference_fasta   = final_refinement.final_assembly_fasta,
            out_basename      = sample_name
  }

  if(call_isnvs) {
    call intrahost.isnvs_per_sample {
        input:
            assembly_fasta = final_refinement.final_assembly_fasta,
            mapped_bam     = final_refinement.aligned_bam
    }
  }

  output {
    # ultimate outputs
    File         final_assembly_fasta        = final_refinement.final_assembly_fasta
    File         aligned_only_reads_bam      = final_refinement.aligned_only_reads_bam
    File         coverage_plot               = final_refinement.coverage_plot
    Int          assembly_length             = final_refinement.assembly_length
    Int          assembly_length_unambiguous = final_refinement.assembly_length_unambiguous
    Int          reads_aligned               = final_refinement.reads_aligned
    Float        mean_coverage               = final_refinement.mean_coverage

    # intermediate outputs
    Array[File?] cleaned_bam               = deplete_taxa.cleaned_bam
    Array[File?] cleaned_fastqc            = deplete_taxa.cleaned_fastqc
    Array[Int?]  depletion_read_count_pre  = deplete_taxa.depletion_read_count_pre
    Array[Int?]  depletion_read_count_post = deplete_taxa.depletion_read_count_post

    Array[File?] taxfilt_bam            = filter_to_taxon.taxfilt_bam
    Array[File?] taxfilt_fastqc         = filter_to_taxon.taxfilt_fastqc
    Array[Int?]  filter_read_count_post = filter_to_taxon.filter_read_count_post

    Array[File]  dedup_bam             = rmdup_ubam.dedup_bam
    Array[File]  dedup_fastqc          = rmdup_ubam.dedup_fastqc
    Array[Int]   dedup_read_count_post = rmdup_ubam.dedup_read_count_post

    File   merged_cleaned_bam = merge_depleted_taxfilt_rmdup.out_bam

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

    File?  isnvsFile = isnvs_per_sample.isnvsFile

    Array[File] aligned_bam                = refine_2x_and_plot.aligned_bam
    Array[File] aligned_only_reads_bam_idx = refine_2x_and_plot.aligned_only_reads_bam_idx
    Array[File] aligned_only_reads_fastqc  = refine_2x_and_plot.aligned_only_reads_fastqc
    Array[File] coverage_tsv               = refine_2x_and_plot.coverage_tsv
    Array[Int]  read_pairs_aligned         = refine_2x_and_plot.read_pairs_aligned
    Array[Float] bases_aligned             = refine_2x_and_plot.bases_aligned

    File merged_refine_input_bam    = merge_refined.out_bam

    File aligned_bam                = final_refinement.aligned_bam
    File aligned_only_reads_bam_idx = final_refinement.aligned_only_reads_bam_idx
    File aligned_only_reads_fastqc  = final_refinement.aligned_only_reads_fastqc
    File coverage_tsv               = final_refinement.coverage_tsv
    Int  read_pairs_aligned         = final_refinement.read_pairs_aligned
    Float bases_aligned             = final_refinement.bases_aligned

    Int    replicate_concordant_sites  = run_discordance.concordant_sites
    Int    replicate_discordant_snps   = run_discordance.discordant_snps
    Int    replicate_discordant_indels = run_discordance.discordant_indels
    Int    num_read_groups             = run_discordance.num_read_groups
    Int    num_libraries               = run_discordance.num_libraries
    File   replicate_discordant_vcf    = run_discordance.discordant_sites_vcf

    Array[String?] deplete_viral_classify_version  = deplete_taxa.viralngs_version
    Array[String?] taxfilt_viral_classify_version  = filter_to_taxon.viralngs_version
    String         assemble_viral_assemble_version = assemble.viralngs_version
    String         scaffold_viral_assemble_version = scaffold.viralngs_version
    Array[String]  refine_viral_assemble_version   = refine_2x_and_plot.viralngs_version
    String         refine_viral_assemble_version   = final_refinement.viralngs_version
    String?        isnvs_viral_phylo_version       = isnvs_per_sample.viralngs_version
  }
}
