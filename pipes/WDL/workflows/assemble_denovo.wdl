version 1.0

import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_intrahost.wdl" as intrahost
import "assemble_refbased.wdl" as assemble_refbased

workflow assemble_denovo {

  meta {
      description: "Assisted de novo viral genome assembly from raw reads."
      author: "Broad Viral Genomics"
      email:  "viral-ngs@broadinstitute.org"
      allowNestedInputs: true
  }

  input {
    File         reads_unmapped_bam

    Array[File]+ reference_genome_fasta

    Array[File]  deplete_bmtaggerDbs = []
    Array[File]  deplete_blastDbs = []
    Array[File]  deplete_bwaDbs =[]

    File?        filter_to_taxon_db
    File         trim_clip_db

    String       sample_name = basename(basename(reads_unmapped_bam, ".bam"), ".cleaned")
  }

  parameter_meta {
    raw_reads_unmapped_bam: { description: "unaligned reads in BAM format", patterns: ["*.bam"] }
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
        lastal_db_fasta    = select_first([filter_to_taxon_db])
    }
  }

  call read_utils.rmdup_ubam {
    input:
      reads_unmapped_bam = select_first([filter_to_taxon.taxfilt_bam, deplete_taxa.cleaned_bam, reads_unmapped_bam])
  }

  call assembly.assemble {
    input:
      reads_unmapped_bam = rmdup_ubam.dedup_bam,
      trim_clip_db       = trim_clip_db,
      always_succeed     = true,
      sample_name        = sample_name
  }

  call assembly.scaffold {
    input:
      contigs_fasta           = assemble.contigs_fasta,
      reads_bam               = select_first([filter_to_taxon.taxfilt_bam, deplete_taxa.cleaned_bam, reads_unmapped_bam]),
      reference_genome_fasta  = reference_genome_fasta
  }

  call assemble_refbased.assemble_refbased {
      input:
          reads_unmapped_bams = [rmdup_ubam.dedup_bam],
          reference_fasta     = scaffold.scaffold_fasta,
          sample_name         = sample_name
  }

  output {
    File    final_assembly_fasta                  = assemble_refbased.assembly_fasta
    File    aligned_only_reads_bam                = assemble_refbased.align_to_self_merged_aligned_only_bam
    File    coverage_plot                         = assemble_refbased.align_to_self_merged_coverage_plot
    Int     assembly_length                       = assemble_refbased.assembly_length
    Int     assembly_length_unambiguous           = assemble_refbased.assembly_length_unambiguous
    Int     reads_aligned                         = assemble_refbased.align_to_self_merged_reads_aligned
    Float   mean_coverage                         = assemble_refbased.align_to_self_merged_mean_coverage
    
    File    cleaned_bam                           = select_first([deplete_taxa.cleaned_bam, reads_unmapped_bam])
    File?   cleaned_fastqc                        = deplete_taxa.cleaned_fastqc
    Int?    depletion_read_count_pre              = deplete_taxa.depletion_read_count_pre
    Int?    depletion_read_count_post             = deplete_taxa.depletion_read_count_post
    
    File?   taxfilt_bam                           = filter_to_taxon.taxfilt_bam
    File?   taxfilt_fastqc                        = filter_to_taxon.taxfilt_fastqc
    Int?    filter_read_count_post                = filter_to_taxon.filter_read_count_post
    
    File    dedup_bam                             = rmdup_ubam.dedup_bam
    File    dedup_fastqc                          = rmdup_ubam.dedup_fastqc
    Int     dedup_read_count_post                 = rmdup_ubam.dedup_read_count_post
    
    File    contigs_fasta                         = assemble.contigs_fasta
    File    subsampBam                            = assemble.subsampBam
    Int     subsample_read_count                  = assemble.subsample_read_count
    
    File    scaffold_fasta                        = scaffold.scaffold_fasta
    File    intermediate_scaffold_fasta           = scaffold.intermediate_scaffold_fasta
    File    intermediate_gapfill_fasta            = scaffold.intermediate_gapfill_fasta
    Int     assembly_preimpute_length             = scaffold.assembly_preimpute_length
    Int     assembly_preimpute_length_unambiguous = scaffold.assembly_preimpute_length_unambiguous
    String  scaffolding_chosen_ref_name           = scaffold.scaffolding_chosen_ref_name
    File    scaffolding_stats                     = scaffold.scaffolding_stats
    File    scaffolding_alt_contigs               = scaffold.scaffolding_alt_contigs
    
    File    isnvsFile                             = assemble_refbased.align_to_self_isnvs_vcf
    
    File    aligned_bam                           = assemble_refbased.aligned_bam
    File    aligned_only_reads_fastqc             = assemble_refbased.align_to_ref_per_input_fastqc[0]
    File    coverage_tsv                          = assemble_refbased.align_to_self_merged_coverage_tsv
    Int     read_pairs_aligned                    = assemble_refbased.align_to_self_merged_read_pairs_aligned
    Float   bases_aligned                         = assemble_refbased.align_to_self_merged_bases_aligned
    
    String? deplete_viral_classify_version        = deplete_taxa.viralngs_version
    String? taxfilt_viral_classify_version        = filter_to_taxon.viralngs_version
    String  assemble_viral_assemble_version       = assemble.viralngs_version
    String  scaffold_viral_assemble_version       = scaffold.viralngs_version
  }

}
