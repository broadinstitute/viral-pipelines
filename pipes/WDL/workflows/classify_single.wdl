version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_utils.wdl" as utils

workflow classify_single {
    meta {
         description: "Runs raw reads through taxonomic classification (Kraken2), human read depletion (based on Kraken2), de novo assembly (SPAdes), and FASTQC/multiQC of reads."
         author: "Broad Viral Genomics"
         email:  "viral-ngs@broadinstitute.org"
         allowNestedInputs: true
    }

    input {
        Array[File]+ reads_bams

        File ncbi_taxdump_tgz

        File spikein_db
        File trim_clip_db

        File kraken2_db_tgz
        File krona_taxonomy_db_kraken2_tgz

        Array[String] taxa_to_dehost         = ["Vertebrata"]
        Array[String] taxa_to_avoid_assembly = ["Vertebrata", "other sequences", "Bacteria"]

        File?         taxid_to_ref_accessions_tsv
    }

    parameter_meta {
        reads_bams: {
          description: "Reads to classify. May be unmapped or mapped or both, paired-end or single-end. Multiple input files will be merged first.",
          patterns: ["*.bam"]
        }
        spikein_db: {
          description: "ERCC spike-in sequences",
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
        krona_taxonomy_db_kraken2_tgz: {
          description: "Krona taxonomy database containing a single file: taxonomy.tab, or possibly just a compressed taxonomy.tab",
          patterns: ["*.tab.zst", "*.tab.gz", "*.tab", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        ncbi_taxdump_tgz: {
          description: "An NCBI taxdump.tar.gz file that contains, at the minimum, a nodes.dmp and names.dmp file.",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        cleaned_fastqc: { 
          description: "Output cleaned fastqc reports in HTML.",
          category: "other"
        }
        deduplicated_reads_unaligned: {
          description: "Deduplication on unaligned reads in BAM format using mvicuna or cdhit.",
          category: "other"
        }
        kraken2_krona_plot: {
          description:"Visualize the results of the Kraken2 analysis with Krona, which disaplys taxonmic hierarchiral data in multi-layerd pie.",
          category:"other"
        }
        kraken2_summary_report: {
          description: "Kraken report output file.",
          category: "other"
        }
        raw_fastqc:{
          description: "Merged raw fastqc reads.",
          category: "other"
        }

    }

    call read_utils.merge_and_reheader_bams as merge_raw_reads {
        input:
            in_bams      = reads_bams
    }
    File reads_bam = merge_raw_reads.out_bam

    call reports.align_and_count as spikein {
        input:
            reads_bam = reads_bam,
            ref_db    = spikein_db
    }
    call metagenomics.kraken2 as kraken2 {
        input:
            reads_bam             = reads_bam,
            kraken2_db_tgz        = kraken2_db_tgz,
            krona_taxonomy_db_tgz = krona_taxonomy_db_kraken2_tgz
    }
    call metagenomics.filter_bam_to_taxa as deplete {
        input:
            classified_bam          = reads_bam,
            classified_reads_txt_gz = kraken2.kraken2_reads_report,
            ncbi_taxonomy_db_tgz    = ncbi_taxdump_tgz,
            exclude_taxa            = true,
            taxonomic_names         = taxa_to_dehost,
            out_filename_suffix     = "hs_depleted"
    }
    call reports.fastqc as fastqc_cleaned {
        input: reads_bam = deplete.bam_filtered_to_taxa
    }
    call metagenomics.filter_bam_to_taxa as filter_acellular {
        input:
            classified_bam          = reads_bam,
            classified_reads_txt_gz = kraken2.kraken2_reads_report,
            ncbi_taxonomy_db_tgz    = ncbi_taxdump_tgz,
            exclude_taxa            = true,
            taxonomic_names         = taxa_to_avoid_assembly,
            out_filename_suffix     = "acellular"
    }
    call read_utils.rmdup_ubam {
       input:
            reads_unmapped_bam = filter_acellular.bam_filtered_to_taxa
    }
    call assembly.assemble as spades {
        input:
            reads_unmapped_bam = rmdup_ubam.dedup_bam,
            trim_clip_db       = trim_clip_db,
            always_succeed     = true
    }
    call metagenomics.report_primary_kraken_taxa {
        input:
            kraken_summary_report = kraken2.kraken2_summary_report
    }

    if(defined(taxid_to_ref_accessions_tsv)) {
      # download (multi-segment) genomes for each reference, fasta filename = colon-concatenated accession list
      scatter(taxon in read_tsv(select_first([taxid_to_ref_accessions_tsv]))) {
          # taxon = [taxid, isolate_prefix, taxname, semicolon_delim_accession_list]
          call utils.string_split {
              input:
                  joined_string = taxon[3],
                  delimiter = ":"
          }
          call ncbi.download_annotations {
              input:
                  accessions = string_split.tokens,
                  combined_out_prefix = sub(taxon[3], ":", "-")  # singularity does not like colons in filenames
          }
      }

      # subset reference genomes to those with ANI hits to contigs and cluster reference hits by any ANI similarity to each other
      call assembly.select_references {
          input:
              reference_genomes_fastas = download_annotations.combined_fasta,
              contigs_fasta = spades.contigs_fasta
      }

      # get taxid and taxname from taxid_to_ref_accessions_tsv
      scatter(top_match in select_references.top_matches_per_cluster_basenames) {
        call utils.fetch_row_from_tsv as tax_lookup {
          input:
              tsv = select_first([taxid_to_ref_accessions_tsv]),
              idx_col = "accessions",
              idx_val = sub(top_match, "-", ":"),
              add_header = ["taxid", "isolate_prefix", "taxname", "accessions"]
        }
        Int skani_hit_taxid = tax_lookup.map["taxid"]
        String skani_hit_taxname = tax_lookup.map["taxname"]
      }
    }


    output {
        File   cleaned_reads_unaligned_bam     = deplete.bam_filtered_to_taxa
        File   deduplicated_reads_unaligned    = rmdup_ubam.dedup_bam
        File   contigs_fasta                   = spades.contigs_fasta
        
        Int    read_counts_raw                 = deplete.classified_taxonomic_filter_read_count_pre
        Int    read_counts_depleted            = deplete.classified_taxonomic_filter_read_count_post
        Int    read_counts_dedup               = rmdup_ubam.dedup_read_count_post
        Int    read_counts_prespades_subsample = spades.subsample_read_count
        
        File   kraken2_summary_report          = kraken2.kraken2_summary_report
        File   kraken2_krona_plot              = kraken2.krona_report_html
        File   kraken2_top_taxa_report         = report_primary_kraken_taxa.ranked_focal_report
        String kraken2_focal_taxon_name        = report_primary_kraken_taxa.focal_tax_name
        Int    kraken2_focal_total_reads       = report_primary_kraken_taxa.total_focal_reads
        String kraken2_top_taxon_id            = report_primary_kraken_taxa.tax_id
        String kraken2_top_taxon_name          = report_primary_kraken_taxa.tax_name
        String kraken2_top_taxon_rank          = report_primary_kraken_taxa.tax_rank
        Int    kraken2_top_taxon_num_reads     = report_primary_kraken_taxa.num_reads
        Float  kraken2_top_taxon_pct_of_focal  = report_primary_kraken_taxa.percent_of_focal

        Int    skani_num_hits                  = length(select_first([select_references.top_matches_per_cluster_basenames, []]))
        File?  skani_contigs_to_refs_dist_tsv  = select_references.skani_dist_full_tsv
        Array[Int]? skani_hits_taxids          = skani_hit_taxid
        Array[String]? skani_hits_taxnames     = skani_hit_taxname

        File   raw_fastqc                      = merge_raw_reads.fastqc
        File   cleaned_fastqc                  = fastqc_cleaned.fastqc_html
        File   spikein_report                  = spikein.report
        String spikein_tophit                  = spikein.top_hit_id
        String spikein_pct_of_total_reads      = spikein.pct_total_reads_mapped
        String spikein_pct_lesser_hits         = spikein.pct_lesser_hits_of_mapped
        
        String kraken2_viral_classify_version  = kraken2.viralngs_version
        String deplete_viral_classify_version  = deplete.viralngs_version
        String spades_viral_assemble_version   = spades.viralngs_version
    }
}
