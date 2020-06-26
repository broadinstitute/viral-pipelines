version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_reports.wdl" as reports

workflow classify_single {
    meta {
         description: "Runs raw reads through taxonomic classification (Kraken2), human read depletion (based on Kraken2), de novo assembly (SPAdes), and FASTQC/multiQC of reads."
         author: "Broad Viral Genomics"
         email:  "viral-ngs@broadinstitute.org"
    }

    input {
        File  reads_bam

        File  ncbi_taxdump_tgz

        File  spikein_db
        File  trim_clip_db

        File  kraken2_db_tgz
        File  krona_taxonomy_db_kraken2_tgz
    }

    parameter_meta {
        reads_bam: {
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
        ncbi_taxdump_tgz: {
          description: "An NCBI taxdump.tar.gz file that contains, at the minimum, a nodes.dmp and names.dmp file.",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
    }

    call reports.fastqc as fastqc_raw {
        input: reads_bam = reads_bam
    }
    call reports.align_and_count as spikein {
        input:
            reads_bam = reads_bam,
            ref_db = spikein_db
    }
    call metagenomics.kraken2 as kraken2 {
        input:
            reads_bam = reads_bam,
            kraken2_db_tgz = kraken2_db_tgz,
            krona_taxonomy_db_tgz = krona_taxonomy_db_kraken2_tgz
    }
    call metagenomics.filter_bam_to_taxa as deplete {
        input:
            classified_bam = reads_bam,
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
            classified_bam = reads_bam,
            classified_reads_txt_gz = kraken2.kraken2_reads_report,
            ncbi_taxonomy_db_tgz = ncbi_taxdump_tgz,
            exclude_taxa = true,
            taxonomic_names = ["Vertebrata", "other sequences", "Bacteria"],
            out_filename_suffix = "acellular"
    }
    call read_utils.rmdup_ubam {
       input:
            reads_unmapped_bam = filter_acellular.bam_filtered_to_taxa
    }
    call assembly.assemble as spades {
        input:
            assembler = "spades",
            reads_unmapped_bam = rmdup_ubam.dedup_bam,
            trim_clip_db = trim_clip_db,
            always_succeed = true
    }

    output {
        File cleaned_reads_unaligned_bam  = deplete.bam_filtered_to_taxa
        File deduplicated_reads_unaligned = rmdup_ubam.dedup_bam
        File contigs_fasta                = spades.contigs_fasta

        Int  read_counts_raw                 = deplete.classified_taxonomic_filter_read_count_pre
        Int  read_counts_depleted            = deplete.classified_taxonomic_filter_read_count_post
        Int  read_counts_dedup               = rmdup_ubam.dedup_read_count_post
        Int  read_counts_prespades_subsample = spades.subsample_read_count

        File kraken2_summary_report = kraken2.kraken2_summary_report
        File kraken2_krona_plot     = kraken2.krona_report_html

        String kraken2_viral_classify_version = kraken2.viralngs_version
        String deplete_viral_classify_version = deplete.viralngs_version
        String spades_viral_assemble_version  = spades.viralngs_version
    }
}
