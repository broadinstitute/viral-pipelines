version 1.0

import "../tasks/tasks_demux.wdl" as demux
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_ncbi.wdl" as ncbi

workflow demux_deplete {
    meta {
        description: "Picard-based demultiplexing and basecalling from a tarball of a raw BCL directory, followed by QC metrics, depletion, and SRA submission prep."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        File         flowcell_tgz
        Array[File]+ samplesheets  ## must be in lane order!

        File?        sample_rename_map
        File?        biosample_map

        File         spikein_db
        Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? bwaDbs
    }

    parameter_meta {
        flowcell_tgz: {
            description: "Illumina BCL directory compressed as tarball. Must contain RunInfo.xml, SampleSheet.csv, RTAComplete.txt, and Data/Intensities/BaseCalls/*",
            patterns: ["*.tar.gz", ".tar.zst", ".tar.bz2", ".tar.lz4", ".tgz"]
        }
        samplesheets: {
            description: "Custom formatted 'extended' format tsv samplesheets that will override any SampleSheet.csv in the illumina BCL directory. Must supply one file per lane of the flowcell, and must provide them in lane order. Required tsv column headings are: sample, library_id_per_sample, barcode_1, barcode_2 (if paired reads, omit if single-end), library_strategy, library_source, library_selection, design_description. 'sample' must correspond to a biological sample. 'sample' x 'library_id_per_sample' must be unique within a samplesheet and correspond to independent libraries from the same original sample. barcode_1 and barcode_2 must correspond to the actual index sequence. Remaining columns must follow strict ontology: see 3rd tab of https://www.ncbi.nlm.nih.gov/core/assets/sra/files/SRA_metadata_acc_example.xlsx for controlled vocabulary and term definitions.",
            patterns: ["*.txt", "*.tsv"]
        }
        sample_rename_map: {
            description: "If 'samples' need to be renamed, provide a two-column tsv that contains at least the following columns: internal_id, external_id. All samples will be renamed prior to analysis. Any samples described in the samplesheets that are not present in sample_rename_map will be unaltered. If this is omitted, no samples will be renamed.",
            patterns: ["*.txt", "*.tsv"]
        }
        biosample_map: {
            description: "A two-column tsv that contains at least the following columns: sample_name, accession. sample_name refers to the external sample id, accession is the NCBI BioSample accession number (SAMNxxx). If this file is omitted, SRA submission prep will be skipped.",
            patterns: ["*.txt", "*.tsv"]
        }
    }

    #### demux each lane (rename samples if requested)
    scatter(lane_sheet in zip(range(length(samplesheets)), samplesheets)) {
        call demux.samplesheet_rename_ids {
            input:
                old_sheet = lane_sheet.right,
                rename_map = sample_rename_map
        }
        call demux.illumina_demux as illumina_demux {
            input:
                flowcell_tgz = flowcell_tgz,
                lane = lane_sheet.left + 1,
                samplesheet = samplesheet_rename_ids.new_sheet
        }
    }

    #### human depletion & spike-in counting for all files
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

    #### SRA submission prep
    if(defined(biosample_map)) {
        call ncbi.sra_meta_prep {
            input:
                cleaned_bam_filepaths = deplete.cleaned_bam,
                biosample_map = select_first([biosample_map]),
                library_metadata = samplesheet_rename_ids.new_sheet,
                platform = "ILLUMINA",
                out_name = "sra_metadata-~{basename(flowcell_tgz, '.tar.gz')}.tsv"
        }
    }

    #### summary stats
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

    # TO DO: flag all libraries where highest spike-in is not what was expected in extended samplesheet

    output {
        Array[File] raw_reads_unaligned_bams     = flatten(illumina_demux.raw_reads_unaligned_bams)
        Array[Int]  read_counts_raw = deplete.depletion_read_count_pre

        Array[File] cleaned_reads_unaligned_bams = deplete.cleaned_bam
        Array[Int]  read_counts_depleted = deplete.depletion_read_count_post

        File?       sra_metadata          = sra_meta_prep.sra_metadata

        Array[File] demux_metrics         = illumina_demux.metrics
        Array[File] demux_commonBarcodes  = illumina_demux.commonBarcodes
        Array[File] demux_outlierBarcodes = illumina_demux.outlierBarcodes

        File        multiqc_report_raw     = multiqc_raw.multiqc_report
        File        multiqc_report_cleaned = multiqc_cleaned.multiqc_report
        File        spikein_counts         = spike_summary.count_summary

        String      demux_viral_core_version = illumina_demux.viralngs_version[0]
    }
}
