version 1.0

#DX_SKIP_WORKFLOW

import "../tasks/tasks_demux.wdl" as demux
import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_reports.wdl" as reports

workflow demux_metag {
    input {
        File spikein_db
        File trim_clip_db
        Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? bwaDbs
    }

    call demux.illumina_demux as illumina_demux

    scatter(raw_reads in illumina_demux.raw_reads_unaligned_bams) {
        call reports.spikein_report as spikein {
            input:
                reads_bam = raw_reads,
                spikein_db = spikein_db
        }
        call taxon_filter.deplete_taxa as deplete {
            input:
                raw_reads_unmapped_bam = raw_reads,
                bmtaggerDbs = bmtaggerDbs,
                blastDbs = blastDbs,
                bwaDbs = bwaDbs
        }
        call assembly.assemble as spades {
            input:
                assembler = "spades",
                reads_unmapped_bam = deplete.cleaned_bam,
                trim_clip_db = trim_clip_db,
                always_succeed = true
        }
        call metagenomics.kaiju as kaiju {
          input:
            reads_unmapped_bam = raw_reads
        }
    }

    call reports.MultiQC as multiqc_raw {
        input:
            input_files = illumina_demux.raw_reads_fastqc_zip,
            file_name   = "multiqc-raw.html"
    }

    call reports.MultiQC as multiqc_cleaned {
        input:
            input_files = deplete.cleaned_fastqc_zip,
            file_name   = "multiqc-cleaned.html"
    }

    call metagenomics.krakenuniq as krakenuniq {
        input:
            reads_unmapped_bam = illumina_demux.raw_reads_unaligned_bams
    }

    call reports.spikein_summary as spike_summary {
        input:
            spikein_count_txt = spikein.report
    }

    call reports.aggregate_metagenomics_reports as metag_summary_report {
        input:
            kraken_summary_reports = krakenuniq.krakenuniq_summary_reports
    }

    output {
        Array[File] raw_reads_unaligned_bams     = illumina_demux.raw_reads_unaligned_bams
        Array[File] cleaned_reads_unaligned_bams = deplete.cleaned_bam
        Array[File] contigs_fastas               = spades.contigs_fasta

        Array[Int]  read_counts_raw = deplete.depletion_read_count_pre
        Array[Int]  read_counts_depleted = deplete.depletion_read_count_post
        Array[Int]  read_counts_prespades_subsample = spades.subsample_read_count

        File        demux_metrics            = illumina_demux.metrics
        File        demux_commonBarcodes     = illumina_demux.commonBarcodes
        File        demux_outlierBarcodes    = illumina_demux.outlierBarcodes

        File        multiqc_report_raw     = multiqc_raw.multiqc_report
        File        multiqc_report_cleaned = multiqc_cleaned.multiqc_report
        File        spikein_counts         = spike_summary.spikein_summary
        File        metagenomics_krona     = krakenuniq.krona_report_merged_html
        File        metagenomics_summary   = metag_summary_report.krakenuniq_aggregate_taxlevel_summary

        String      demux_viral_core_version          = illumina_demux.viralngs_version
        String      krakenuniq_viral_classify_version = krakenuniq.viralngs_version
        String      deplete_viral_classify_version    = deplete.viralngs_version[0]
        String      spades_viral_assemble_version     = spades.viralngs_version[0]
    }
}
