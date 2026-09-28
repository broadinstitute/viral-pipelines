version 1.1

import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_reports.wdl" as reports

workflow classify_multi {
    meta {
         description: "Runs raw reads through taxonomic classification (Kraken2), human read depletion (based on Kraken2), de novo assembly (SPAdes), and FASTQC/multiQC of reads."
         author: "Broad Viral Genomics"
         email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+ reads_bams

        File ncbi_taxdump_tgz

        File spikein_db
        File trim_clip_db

        File kraken2_db_tgz
        File krona_taxonomy_db_kraken2_tgz

        Int min_reads_for_rmdup    =  5000000
        Int max_reads_for_assembly = 10000000
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
        ncbi_taxdump_tgz: {
          description: "An NCBI taxdump.tar.gz file that contains, at the minimum, a nodes.dmp and names.dmp file.",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        min_reads_for_rmdup: {
          description: "Read sets smaller than this skip kmer-depth normalization and are passed to assembly unchanged."
        }
        max_reads_for_assembly: {
          description: "Cap on read pairs passed to de novo assembly. If normalization leaves more than this, the read set is randomly downsampled to this count."
        }
    }

    scatter(raw_reads in reads_bams) {
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
                reads_bam             = raw_reads,
                kraken2_db_tgz        = kraken2_db_tgz,
                krona_taxonomy_db_tgz = krona_taxonomy_db_kraken2_tgz
        }
        call metagenomics.filter_bam_to_taxa as deplete {
            input:
                classified_bam          = raw_reads,
                classified_reads_txt_gz = kraken2.kraken2_reads_report,
                ncbi_taxonomy_db_tgz    = ncbi_taxdump_tgz,
                exclude_taxa            = true,
                taxonomic_names         = ["Vertebrata"],
                out_filename_suffix     = "hs_depleted"
        }
        call metagenomics.filter_bam_to_taxa as filter_acellular {
            input:
                classified_bam          = raw_reads,
                classified_reads_txt_gz = kraken2.kraken2_reads_report,
                ncbi_taxonomy_db_tgz    = ncbi_taxdump_tgz,
                exclude_taxa            = true,
                taxonomic_names         = ["Vertebrata", "other sequences", "Bacteria"],
                out_filename_suffix     = "acellular"
        }
    }

    scatter(clean_reads in filter_acellular.bam_filtered_to_taxa) {
        call read_utils.bbnorm_bam {
           input:
                reads_bam        = clean_reads,
                min_input_reads  = min_reads_for_rmdup,
                max_output_reads = max_reads_for_assembly
        }
        call assembly.assemble as spades {
            input:
                reads_unmapped_bam = bbnorm_bam.bbnorm_bam,
                trim_clip_db       = trim_clip_db,
                always_succeed     = true
        }
    }

    call reports.multiqc_from_bams as multiqc_raw {
        input:
            input_bams   = reads_bams,
            out_basename = "multiqc-raw"
    }

    call reports.multiqc_from_bams as multiqc_cleaned {
        input:
            input_bams   = deplete.bam_filtered_to_taxa,
            out_basename = "multiqc-cleaned"
    }

    call reports.multiqc_from_bams as multiqc_assembly_input {
        input:
            input_bams   = bbnorm_bam.bbnorm_bam,
            out_basename = "multiqc-assembly-input"
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
            reports_txt_gz        = kraken2.kraken2_summary_report,
            krona_taxonomy_db_tgz = krona_taxonomy_db_kraken2_tgz,
            input_type            = "kraken2",
            out_basename          = "merged-kraken2.krona"
    }

    output {
        Array[File] cleaned_reads_unaligned_bams    = deplete.bam_filtered_to_taxa
        Array[File] reads_assembly_input_ubams      = bbnorm_bam.bbnorm_bam
        Array[File] contigs_fastas                  = spades.contigs_fasta
        
        Array[Int]  read_counts_raw                 = deplete.classified_taxonomic_filter_read_count_pre
        Array[Int]  read_counts_depleted            = deplete.classified_taxonomic_filter_read_count_post
        Array[Int]  read_counts_assembly_input      = bbnorm_bam.bbnorm_read_count_post
        Array[Int]  read_counts_prespades_subsample = spades.subsample_read_count
        
        File        multiqc_report_raw              = multiqc_raw.multiqc_report
        File        multiqc_report_cleaned          = multiqc_cleaned.multiqc_report
        File        multiqc_report_assembly_input   = multiqc_assembly_input.multiqc_report
        File        spikein_counts                  = spike_summary.count_summary
        File        kraken2_merged_krona            = krona_merge_kraken2.krona_report_html
        File        kraken2_summary                 = metag_summary_report.krakenuniq_aggregate_taxlevel_summary
        
        Array[File] kraken2_summary_reports         = kraken2.kraken2_summary_report
        Array[File] kraken2_krona_by_sample         = kraken2.krona_report_html
        
        String      kraken2_viral_classify_version  = kraken2.viralngs_version[0]
        String      deplete_viral_classify_version  = deplete.viralngs_version[0]
        String      spades_viral_assemble_version   = spades.viralngs_version[0]
    }
}
