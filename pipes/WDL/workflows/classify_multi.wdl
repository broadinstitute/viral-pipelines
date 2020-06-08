version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_reports.wdl" as reports

workflow classify_multi {
    meta {
         description: "Runs raw reads through taxonomic classification (Kraken2), human read depletion (based on Kraken2), de novo assembly (SPAdes), taxonomic classification of contigs (BLASTx), and FASTQC/multiQC of reads."
         author: "Broad Viral Genomics"
         email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+ reads_bams

        File  ncbi_taxdump_tgz

        File  spikein_db
        File  trim_clip_db

        File  kraken2_db_tgz
        File  krona_taxonomy_db_kraken2_tgz
        File? blast_db_tgz
        File? krona_taxonomy_db_blast_tgz
    }

    parameter_meta {
        reads_bams: {
          description: "Reads to classify. May be unmapped or mapped or both, paired-end or single-end.",
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
        blast_db_tgz: {
          description: "Pre-built BLAST database tarball containing an indexed blast database named 'nr'",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        krona_taxonomy_db_blast_tgz: {
          description: "Krona taxonomy database: a tarball containing a taxonomy.tab file as well as accession to taxid mapping (a kraken-based taxonomy database will not suffice).",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        ncbi_taxdump_tgz: {
          description: "An NCBI taxdump.tar.gz file that contains, at the minimum, a nodes.dmp and names.dmp file.",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
    }

    scatter(raw_reads in reads_bams) {
        call reports.fastqc as fastqc_raw {
            input: reads_bam = raw_reads
        }
        call reports.align_and_count as spikein {
            input:
                reads_bam = raw_reads,
                ref_db = spikein_db
        }
    }

    scatter(raw_reads in reads_bams) {
        # separate scatter blocks speeds up the gathers in DNAnexus and provides independent failure blocks
        call metagenomics.kraken2 as kraken2 {
            input:
                reads_bam = raw_reads,
                kraken2_db_tgz = kraken2_db_tgz,
                krona_taxonomy_db_tgz = krona_taxonomy_db_kraken2_tgz
        }
        call metagenomics.filter_bam_to_taxa as deplete {
            input:
                classified_bam = raw_reads,
                classified_reads_txt_gz = kraken2.kraken2_reads_report,
                ncbi_taxonomy_db_tgz = ncbi_taxdump_tgz,
                exclude_taxa = true,
                taxonomic_names = ["Vertebrata"],
                out_filename_suffix = "hs_depleted"
        }
        call reports.fastqc as fastqc_cleaned {
            input: reads_bam = deplete.bam_filtered_to_taxa
        }
        call metagenomics.filter_bam_to_taxa as filter_acellular {
            input:
                classified_bam = raw_reads,
                classified_reads_txt_gz = kraken2.kraken2_reads_report,
                ncbi_taxonomy_db_tgz = ncbi_taxdump_tgz,
                exclude_taxa = true,
                taxonomic_names = ["Vertebrata", "other sequences", "Bacteria"],
                out_filename_suffix = "acellular"
        }
    }

    scatter(clean_reads in filter_acellular.bam_filtered_to_taxa) {
        call read_utils.rmdup_ubam {
           input:
                reads_unmapped_bam = clean_reads
        }
        call assembly.assemble as spades {
            input:
                assembler = "spades",
                reads_unmapped_bam = rmdup_ubam.dedup_bam,
                trim_clip_db = trim_clip_db,
                spades_min_contig_len = 800,
                always_succeed = true
        }
        if(defined(blast_db_tgz) && defined(krona_taxonomy_db_blast_tgz)) {
            call metagenomics.blastx as blastx {
                input:
                    contigs_fasta = spades.contigs_fasta,
                    blast_db_tgz = select_first([blast_db_tgz]),
                    krona_taxonomy_db_tgz = select_first([krona_taxonomy_db_blast_tgz])
            }
        }
    }

    call reports.MultiQC as multiqc_raw {
        input:
            input_files = fastqc_raw.fastqc_zip,
            file_name   = "multiqc-raw.html"
    }

    call reports.MultiQC as multiqc_cleaned {
        input:
            input_files = fastqc_cleaned.fastqc_zip,
            file_name   = "multiqc-cleaned.html"
    }

    call reports.MultiQC as multiqc_dedup {
        input:
            input_files = rmdup_ubam.dedup_fastqc_zip,
            file_name   = "multiqc-dedup.html"
    }

    call reports.align_and_count_summary as spike_summary {
        input:
            counts_txt = spikein.report
    }

    call reports.aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = kraken2.kraken2_summary_report
    }

    call metagenomics.krona as krona_merge_kraken2 {
        input:
            reports_txt_gz = kraken2.kraken2_summary_report,
            krona_taxonomy_db_tgz = krona_taxonomy_db_kraken2_tgz,
            input_type = "kraken2",
            out_basename = "merged-kraken2.krona"
    }

    if(defined(blast_db_tgz) && defined(krona_taxonomy_db_blast_tgz)) {
        call metagenomics.krona_merge as krona_merge_blastx {
            input:
                krona_reports = select_all(blastx.krona_report_html),
                out_basename = "merged-spades-blastx.krona"
        }
    }

    output {
        Array[File] cleaned_reads_unaligned_bams = deplete.bam_filtered_to_taxa
        Array[File] deduplicated_reads_unaligned = rmdup_ubam.dedup_bam
        Array[File] contigs_fastas               = spades.contigs_fasta

        Array[Int]  read_counts_raw                 = deplete.classified_taxonomic_filter_read_count_pre
        Array[Int]  read_counts_depleted            = deplete.classified_taxonomic_filter_read_count_post
        Array[Int]  read_counts_dedup               = rmdup_ubam.dedup_read_count_post
        Array[Int]  read_counts_prespades_subsample = spades.subsample_read_count

        File        multiqc_report_raw     = multiqc_raw.multiqc_report
        File        multiqc_report_cleaned = multiqc_cleaned.multiqc_report
        File        multiqc_report_dedup   = multiqc_dedup.multiqc_report
        File        spikein_counts         = spike_summary.count_summary
        File        kraken2_merged_krona   = krona_merge_kraken2.krona_report_html
        File        kraken2_summary        = metag_summary_report.krakenuniq_aggregate_taxlevel_summary
        File?       blastx_merged_krona   = krona_merge_blastx.krona_report_html

        Array[File] kraken2_summary_reports = kraken2.kraken2_summary_report
        Array[File] kraken2_krona_by_sample = kraken2.krona_report_html
        Array[File] blastx_report_by_sample = select_all(blastx.blast_report)
        Array[File] blastx_krona_by_sample  = select_all(blastx.krona_report_html)

        String      kraken2_viral_classify_version = kraken2.viralngs_version[0]
        String      deplete_viral_classify_version    = deplete.viralngs_version[0]
        String      spades_viral_assemble_version     = spades.viralngs_version[0]
    }
}
