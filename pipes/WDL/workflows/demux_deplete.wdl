version 1.0

import "../tasks/tasks_demux.wdl" as demux
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_reports.wdl" as reports

workflow demux_plus {
    meta {
        description: "Picard-based demultiplexing and basecalling from a tarball of a raw BCL directory, followed by QC metrics and depletion. Intended for automatic triggering post upload on DNAnexus."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        File         flowcell_tgz
        Array[File]+ samplesheets

        File spikein_db
        Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? bwaDbs
    }

    scatter(lane_sheet in zip(range(length(samplesheets)), samplesheets)) {
        call demux.illumina_demux as illumina_demux {
            input:
                flowcell_tgz = flowcell_tgz,
                lane = lane_sheet.left,
                samplesheet = lane_sheet.right
        }

    }

    scatter(raw_reads in flatten(illumina_demux.raw_reads_unaligned_bams)) {
        call reports.align_and_count as spikein {
            input:
                reads_bam = raw_reads,
                ref_db = spikein_db
        }
        call taxon_filter.deplete_taxa as deplete {
            input:
                raw_reads_unmapped_bam = raw_reads,
                bmtaggerDbs = bmtaggerDbs,
                blastDbs = blastDbs,
                bwaDbs = bwaDbs
        }
    }

    call reports.MultiQC as multiqc_raw {
        input:
            input_files = flatten(illumina_demux.raw_reads_fastqc_zip),
            file_name   = "multiqc-raw.html"
    }

    call reports.MultiQC as multiqc_cleaned {
        input:
            input_files = deplete.cleaned_fastqc_zip,
            file_name   = "multiqc-cleaned.html"
    }

    call reports.align_and_count_summary as spike_summary {
        input:
            counts_txt = spikein.report
    }

    output {
        Array[File] raw_reads_unaligned_bams     = flatten(illumina_demux.raw_reads_unaligned_bams)
        Array[File] cleaned_reads_unaligned_bams = deplete.cleaned_bam

        Array[Int]  read_counts_raw = deplete.depletion_read_count_pre
        Array[Int]  read_counts_depleted = deplete.depletion_read_count_post

        Array[File] demux_metrics            = illumina_demux.metrics
        Array[File] demux_commonBarcodes     = illumina_demux.commonBarcodes
        Array[File] demux_outlierBarcodes    = illumina_demux.outlierBarcodes

        File        multiqc_report_raw     = multiqc_raw.multiqc_report
        File        multiqc_report_cleaned = multiqc_cleaned.multiqc_report
        File        spikein_counts         = spike_summary.count_summary

        String      demux_viral_core_version          = illumina_demux.viralngs_version[0]
    }
}
