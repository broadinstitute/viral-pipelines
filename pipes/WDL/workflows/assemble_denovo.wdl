version 1.0

import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_ncbi.wdl" as ncbi
import "assemble_refbased.wdl" as assemble_refbased

workflow assemble_denovo {

  meta {
      description: "Assisted de novo viral genome assembly from raw reads."
      author: "Broad Viral Genomics"
      email:  "viral-ngs@broadinstitute.org"
      allowNestedInputs: true
  }

  input {
    Array[File]+ reads_unmapped_bams

    Array[File]+ reference_genome_fasta

    Array[File]  deplete_bmtaggerDbs = []
    Array[File]  deplete_blastDbs = []
    Array[File]  deplete_bwaDbs =[]

    File?        filter_to_taxon_db
    File         trim_clip_db

    String       out_basename = basename(basename(reads_unmapped_bams[0], ".bam"), ".cleaned")
    String?      sample_original_name
  }

  parameter_meta {
    raw_reads_unmapped_bams: { description: "unaligned reads in BAM format", patterns: ["*.bam"] }
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
    out_basename: { description: "a filename-friendly basename for output files" }
    sample_original_name: { description: "a (possibly filename-unfriendly) sample name for fasta and bam headers" }
  }

  # parallelize across provided input read files
  scatter(reads_unmapped_bam in reads_unmapped_bams) {

    # rename SM value in bam header if requested
    if(defined(sample_original_name)) {
      call read_utils.merge_and_reheader_bams as renamed_reads {
          input:
              in_bams      = [reads_unmapped_bam],
              sample_name  = sample_original_name,
              out_basename = out_basename
      }
    }
    File reads_unmapped_renamed_bams = select_first([renamed_reads.out_bam, reads_unmapped_bam])

    # deplete host if requested
    if(length(deplete_bmtaggerDbs) + length(deplete_blastDbs) + length(deplete_bwaDbs) > 0) {
      call taxon_filter.deplete_taxa {
        input:
          raw_reads_unmapped_bam = reads_unmapped_renamed_bams,
          bmtaggerDbs            = deplete_bmtaggerDbs,
          blastDbs               = deplete_blastDbs,
          bwaDbs                 = deplete_bwaDbs
      }
    }
    File reads_depleted_bams = select_first([deplete_taxa.cleaned_bam, reads_unmapped_bam])

    # select reads if requested
    if(defined(filter_to_taxon_db)) {
      call taxon_filter.filter_to_taxon {
        input:
          reads_unmapped_bam = reads_depleted_bams,
          lastal_db_fasta    = select_first([filter_to_taxon_db])
      }
    }
    File reads_taxfilt_bams = select_first([filter_to_taxon.taxfilt_bam, reads_depleted_bams])

    # alignment-free PCR duplicate removal
    call read_utils.rmdup_ubam {
      input:
        reads_unmapped_bam = reads_taxfilt_bams
    }
  }

  # merge all reads into single file
  call read_utils.merge_and_reheader_bams as merge_dedup_reads {
      input:
          in_bams      = rmdup_ubam.dedup_bam,
          out_basename = out_basename
  }
  call read_utils.merge_and_reheader_bams as merge_cleaned_reads {
      input:
          in_bams      = reads_depleted_bams,
          out_basename = out_basename
  }
  call read_utils.merge_and_reheader_bams as merge_taxfilt_reads {
      input:
          in_bams      = reads_taxfilt_bams,
          out_basename = out_basename
  }

  # denovo assembly pipeline below
  call assembly.assemble {
    input:
      reads_unmapped_bam = merge_dedup_reads.out_bam,
      trim_clip_db       = trim_clip_db,
      always_succeed     = true,
      sample_name        = out_basename
  }

  call assembly.scaffold {
    input:
      contigs_fasta           = assemble.contigs_fasta,
      reads_bam               = merge_dedup_reads.out_bam,
      reference_genome_fasta  = reference_genome_fasta
  }

  call assemble_refbased.assemble_refbased as refine {
      input:
          reads_unmapped_bams = reads_depleted_bams, # assemble_refbased will scatter on individual bams
          reference_fasta     = scaffold.scaffold_fasta,
          sample_name         = out_basename
  }

  call assemble_refbased.assemble_refbased as refine2 {
      input:
          reads_unmapped_bams = reads_depleted_bams, # assemble_refbased will scatter on individual bams
          reference_fasta     = refine.assembly_fasta,
          sample_name         = out_basename
  }

  if (defined(sample_original_name)) {
    call ncbi.rename_fasta_header {
      input:
        genome_fasta = refine2.assembly_fasta,
        new_name     = select_first([sample_original_name])
    }
  }

  output {
    File    final_assembly_fasta                  = select_first([rename_fasta_header.renamed_fasta, refine2.assembly_fasta])
    File    aligned_only_reads_bam                = refine2.align_to_self_merged_aligned_only_bam
    File    coverage_plot                         = refine2.align_to_self_merged_coverage_plot
    Int     assembly_length                       = refine2.assembly_length
    Int     assembly_length_unambiguous           = refine2.assembly_length_unambiguous
    Int     reads_aligned                         = refine2.align_to_self_merged_reads_aligned
    Float   mean_coverage                         = refine2.align_to_self_merged_mean_coverage
    
    File    cleaned_bam                           = merge_cleaned_reads.out_bam
    File    cleaned_fastqc                        = merge_cleaned_reads.fastqc
    Int     depletion_read_count_post             = merge_cleaned_reads.read_count
    
    File    taxfilt_bam                           = merge_taxfilt_reads.out_bam
    File    taxfilt_fastqc                        = merge_taxfilt_reads.fastqc
    Int     filter_read_count_post                = merge_taxfilt_reads.read_count
    
    File    dedup_bam                             = merge_dedup_reads.out_bam
    File    dedup_fastqc                          = merge_dedup_reads.fastqc
    Int     dedup_read_count_post                 = merge_dedup_reads.read_count
    
    File    contigs_fasta                         = assemble.contigs_fasta
    File    subsampBam                            = assemble.subsampBam
    Int     subsample_read_count                  = assemble.subsample_read_count
    
    File    scaffold_fasta                        = scaffold.scaffold_fasta
    File    intermediate_scaffold_fasta           = scaffold.intermediate_scaffold_fasta
    File    intermediate_gapfill_fasta            = scaffold.intermediate_gapfill_fasta
    Int     assembly_preimpute_length             = scaffold.assembly_preimpute_length
    Int     assembly_preimpute_length_unambiguous = scaffold.assembly_preimpute_length_unambiguous
    Array[String]  scaffolding_chosen_ref_names   = scaffold.scaffolding_chosen_ref_names
    File    scaffolding_stats                     = scaffold.scaffolding_stats
    File    scaffolding_alt_contigs               = scaffold.scaffolding_alt_contigs

    Int     replicate_concordant_sites            = refine2.replicate_concordant_sites
    Int     replicate_discordant_snps             = refine2.replicate_discordant_snps
    Int     replicate_discordant_indels           = refine2.replicate_discordant_indels
    Int     num_read_groups                       = refine2.num_read_groups
    Int     num_libraries                         = refine2.num_libraries
    File    replicate_discordant_vcf              = refine2.replicate_discordant_vcf

    File    isnvs_vcf                             = refine2.align_to_self_isnvs_vcf
    
    File    aligned_bam                           = refine2.align_to_self_merged_aligned_only_bam
    File    aligned_only_reads_fastqc             = refine2.align_to_ref_fastqc
    File    coverage_tsv                          = refine2.align_to_self_merged_coverage_tsv
    Int     read_pairs_aligned                    = refine2.align_to_self_merged_read_pairs_aligned
    Float   bases_aligned                         = refine2.align_to_self_merged_bases_aligned
    
    String  assembly_method = "viral-ngs/assemble_denovo"
    String  assemble_viral_assemble_version       = assemble.viralngs_version
    String  scaffold_viral_assemble_version       = scaffold.viralngs_version
  }
}
