version 1.0

import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_reports.wdl" as reports
import "assemble_refbased.wdl" as assemble_refbased

workflow metagenomic_denovo {

  meta {
      description: "Assisted de novo viral genome assembly (SPAdes, scaffolding, and polishing) from metagenomic raw reads. Runs raw reads through taxonomic classification (Kraken2), human read depletion (based on Kraken2 and optionally using BWA, BLASTN, and/or BMTAGGER databases), and FASTQC/multiQC of reads."
      author: "Broad Viral Genomics"
      email:  "viral-ngs@broadinstitute.org"
      allowNestedInputs: true
  }

  input {
    String        sample_name
    File          fastq_1
    File?         fastq_2
    String        library_name="1"
    String?       run_date_iso
    String        sequencing_platform

    Array[File]+  reference_genome_fasta

    File          kraken2_db_tgz
    File          krona_taxonomy_tab
    File          ncbi_taxdump_tgz

    Array[File]   deplete_bmtaggerDbs = []
    Array[File]   deplete_blastDbs = []
    Array[File]   deplete_bwaDbs =[]

    Array[String] taxa_to_dehost = ["Vertebrata"]
    Array[String] taxa_to_avoid_assembly = ["Vertebrata", "other sequences", "Bacteria"]

    File?         filter_to_taxon_db
    File?         spikein_db

    File          trim_clip_db
  }

  parameter_meta {
    fastq_1: { description: "Unaligned read1 file in fastq format", patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"] }
    fastq_2: { description: "Unaligned read2 file in fastq format. This should be empty for single-end read conversion and required for paired-end reads. If provided, it must match fastq_1 in length and order.", patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"] }
    sample_name: { description: "Sample name. This is required and will populate the 'SM' read group value and will be used as the output filename (must be filename-friendly)." }
    sequencing_platform: { description: "Sequencing platform. This is required and will populate the 'PL' read group value. Must be one of CAPILLARY, DNBSEQ, HELICOS, ILLUMINA, IONTORRENT, LS454, ONT, PACBIO, or SOLID." }

    reference_genome_fasta: {
      description: "After denovo assembly, large contigs are scaffolded against a reference genome to determine orientation and to join contigs together, before further polishing by reads. You must supply at least one reference genome (all segments/chromomes in a single fasta file). If more than one reference is provided, contigs will be scaffolded against all of them and the one with the most complete assembly will be chosen for downstream polishing.",
      patterns: ["*.fasta"]
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
    spikein_db: {
      description: "ERCC/SDSI spike-in sequences",
      patterns: ["*.fasta", "*.fasta.gz", "*.fasta.zst"]
    }
    trim_clip_db: {
      description: "Adapter sequences to remove via trimmomatic prior to SPAdes assembly",
      patterns: ["*.fasta", "*.fasta.gz", "*.fasta.zst"]
    }
    kraken2_db_tgz: {
      description: "Pre-built Kraken database tarball containing three files: hash.k2d, opts.k2d, and taxo.k2d.",
      patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    krona_taxonomy_tab: {
      description: "Krona taxonomy database containing a single file: taxonomy.tab, or possibly just a compressed taxonomy.tab",
      patterns: ["*.tab.zst", "*.tab.gz", "*.tab", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    ncbi_taxdump_tgz: {
      description: "An NCBI taxdump.tar.gz file that contains, at the minimum, a nodes.dmp and names.dmp file.",
      patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
  }

  # bundle 1 or 2 input files along with their metadata
  call read_utils.FastqToUBAM {
    input:
      fastq_1 = fastq_1,
      fastq_2 = fastq_2,
      sample_name = sample_name,
      library_name = library_name,
      run_date = run_date_iso,
      platform_name = sequencing_platform
  }
  File reads_bam = FastqToUBAM.unmapped_bam

  # metagenomics, QC, etc
  call reports.fastqc as fastqc_raw {
      input: reads_bam = reads_bam
  }
  if(defined(spikein_db)) {
    call reports.align_and_count as spikein {
        input:
            reads_bam = reads_bam,
            ref_db    = select_first([spikein_db])
    }
  }
  call metagenomics.kraken2 as kraken2 {
      input:
          reads_bam             = reads_bam,
          kraken2_db_tgz        = kraken2_db_tgz,
          krona_taxonomy_db_tgz = krona_taxonomy_tab
  }

  # deplete host reads: kraken2 followed by (optional) bmtagger, blastn, and/or bwa
  # resulting read set can be published to SRA
  call metagenomics.filter_bam_to_taxa as deplete_k2 {
      input:
          classified_bam          = reads_bam,
          classified_reads_txt_gz = kraken2.kraken2_reads_report,
          ncbi_taxonomy_db_tgz    = ncbi_taxdump_tgz,
          exclude_taxa            = true,
          taxonomic_names         = taxa_to_dehost,
          out_filename_suffix     = "host_depleted"
  }
  if(length(deplete_bmtaggerDbs) + length(deplete_blastDbs) + length(deplete_bwaDbs) > 0) {
    call taxon_filter.deplete_taxa {
      input:
        raw_reads_unmapped_bam = deplete_k2.bam_filtered_to_taxa,
        bmtaggerDbs            = deplete_bmtaggerDbs,
        blastDbs               = deplete_blastDbs,
        bwaDbs                 = deplete_bwaDbs
    }
  }
  File dehosted_bam = select_first([deplete_taxa.cleaned_bam, deplete_k2.bam_filtered_to_taxa])
  call reports.fastqc as fastqc_dehosted {
      input: reads_bam = dehosted_bam
  }

  # taxonomically focus/filter read set: remove Bacterial and artificial taxa (kraken2) and optionally filter to desired LASTAL database
  call metagenomics.filter_bam_to_taxa as filter_acellular {
      input:
          classified_bam          = dehosted_bam,
          classified_reads_txt_gz = kraken2.kraken2_reads_report,
          ncbi_taxonomy_db_tgz    = ncbi_taxdump_tgz,
          exclude_taxa            = true,
          taxonomic_names         = taxa_to_avoid_assembly,
          out_filename_suffix     = "acellular"
  }
  if(defined(filter_to_taxon_db)) {
    call taxon_filter.filter_to_taxon {
      input:
        reads_unmapped_bam = dehosted_bam,
        lastal_db_fasta    = select_first([filter_to_taxon_db])
    }
  }
  File taxfiltered_bam = select_first([filter_to_taxon.taxfilt_bam, dehosted_bam])
  call reports.fastqc as fastqc_taxfilt {
      input: reads_bam = taxfiltered_bam
  }

  # alignment-free duplicate removal
  call read_utils.rmdup_ubam {
    input:
      reads_unmapped_bam = taxfiltered_bam
  }

  # denovo assembly (with taxfiltered/rmdup reads)
  call assembly.assemble {
    input:
      reads_unmapped_bam = rmdup_ubam.dedup_bam,
      trim_clip_db       = trim_clip_db,
      always_succeed     = true,
      sample_name        = sample_name
  }

  # scaffold contigs to one or more ref genomes and impute with reference
  call assembly.scaffold {
    input:
      contigs_fasta           = assemble.contigs_fasta,
      reads_bam               = dehosted_bam,
      sample_name             = sample_name,
      reference_genome_fasta  = reference_genome_fasta
  }

  # polish/refine with dehosted reads
  call assemble_refbased.assemble_refbased as refine {
      input:
          reads_unmapped_bams = [dehosted_bam],
          reference_fasta     = scaffold.scaffold_fasta,
          sample_name         = sample_name
  }

  output {
    File    final_assembly_fasta                  = refine.assembly_fasta
    File    aligned_only_reads_bam                = refine.align_to_self_merged_aligned_only_bam
    File    coverage_plot                         = refine.align_to_self_merged_coverage_plot
    Int     assembly_length                       = refine.assembly_length
    Int     assembly_length_unambiguous           = refine.assembly_length_unambiguous
    Int     reads_aligned                         = refine.align_to_self_merged_reads_aligned
    Float   mean_coverage                         = refine.align_to_self_merged_mean_coverage

    File    kraken2_summary_report                = kraken2.kraken2_summary_report
    File    kraken2_krona_plot                    = kraken2.krona_report_html

    File    raw_unmapped_bam                      = reads_bam
    File    depleted_bam                          = dehosted_bam
    File    taxfilt_bam                           = taxfiltered_bam
    File    dedup_bam                             = rmdup_ubam.dedup_bam
    File    denovo_in_bam                         = assemble.subsampBam
        
    Int     read_counts_raw                       = deplete_k2.classified_taxonomic_filter_read_count_pre
    Int     read_counts_depleted                  = select_first([deplete_taxa.depletion_read_count_post, deplete_k2.classified_taxonomic_filter_read_count_post])
    Int     read_counts_taxfilt                   = select_first([filter_to_taxon.filter_read_count_post, filter_acellular.classified_taxonomic_filter_read_count_post])
    Int     read_counts_dedup                     = rmdup_ubam.dedup_read_count_post
    Int     read_counts_denovo_in                 = assemble.subsample_read_count

    File    raw_fastqc                            = fastqc_raw.fastqc_html
    File    depleted_fastqc                       = fastqc_dehosted.fastqc_html
    File    taxfilt_fastqc                        = fastqc_taxfilt.fastqc_html

    File    contigs_fasta                         = assemble.contigs_fasta
    File    scaffold_fasta                        = scaffold.scaffold_fasta
    File    intermediate_scaffold_fasta           = scaffold.intermediate_scaffold_fasta
    File    intermediate_gapfill_fasta            = scaffold.intermediate_gapfill_fasta
    Int     assembly_preimpute_length             = scaffold.assembly_preimpute_length
    Int     assembly_preimpute_length_unambiguous = scaffold.assembly_preimpute_length_unambiguous
    Array[String]  scaffolding_chosen_ref_names   = scaffold.scaffolding_chosen_ref_names
    File    scaffolding_stats                     = scaffold.scaffolding_stats
    File    scaffolding_alt_contigs               = scaffold.scaffolding_alt_contigs

    Int     replicate_concordant_sites            = refine.replicate_concordant_sites
    Int     replicate_discordant_snps             = refine.replicate_discordant_snps
    Int     replicate_discordant_indels           = refine.replicate_discordant_indels
    Int     num_read_groups                       = refine.num_read_groups
    Int     num_libraries                         = refine.num_libraries
    File    replicate_discordant_vcf              = refine.replicate_discordant_vcf

    File    isnvs_vcf                             = refine.align_to_self_isnvs_vcf
    
    File    aligned_bam                           = refine.align_to_self_merged_aligned_only_bam
    File    coverage_tsv                          = refine.align_to_self_merged_coverage_tsv
    Int     read_pairs_aligned                    = refine.align_to_self_merged_read_pairs_aligned
    Float   bases_aligned                         = refine.align_to_self_merged_bases_aligned

    File?   spikein_hits                          = spikein.report
    String? spikein_tophit                        = spikein.top_hit_id
    String? spikein_pct_of_total_reads            = spikein.pct_total_reads_mapped
    String? spikein_pct_lesser_hits               = spikein.pct_lesser_hits_of_mapped

    String  viral_classify_version                = kraken2.viralngs_version
    String  viral_assemble_version                = assemble.viralngs_version
  }
}
