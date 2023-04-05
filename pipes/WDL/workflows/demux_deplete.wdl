version 1.0

import "../tasks/tasks_demux.wdl" as demux
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_ncbi.wdl" as ncbi

workflow demux_deplete {
    meta {
        description: "Picard-based demultiplexing and basecalling from a tarball of a raw BCL directory, followed by QC metrics, depletion, and SRA submission prep."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File         flowcell_tgz
        Array[File]+ samplesheets  ## must be in lane order!
        String?      read_structure

        Boolean      sort_reads=true

        File?        sample_rename_map
        File?        biosample_map
        Int          min_reads_per_bam = 100

        String?      instrument_model_user_specified
        String?      sra_title

        File         spikein_db
        Array[File]? bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]? bwaDbs

        Array[String] default_sample_keys = ["amplicon_set", "control", "batch_lib", "viral_ct"]
        Array[String] default_filename_keys = ["spike_in", "batch_lib"]
    }

    parameter_meta {
        flowcell_tgz: {
            description: "Illumina BCL directory compressed as tarball. Must contain RunInfo.xml, SampleSheet.csv, RTAComplete.txt, and Data/Intensities/BaseCalls/*",
            patterns: ["*.tar.gz", ".tar.zst", ".tar.bz2", ".tar.lz4", ".tgz"]
            category: "required"
            
        }
        samplesheets: {
            description: "Custom formatted 'extended' format tsv samplesheets that will override any SampleSheet.csv in the illumina BCL directory. Must supply one file per lane of the flowcell, and must provide them in lane order. Required tsv column headings are: sample, library_id_per_sample, barcode_1, barcode_2 (if paired reads, omit if single-end), library_strategy, library_source, library_selection, design_description. 'sample' must correspond to a biological sample. 'sample' x 'library_id_per_sample' must be unique within a samplesheet and correspond to independent libraries from the same original sample. barcode_1 and barcode_2 must correspond to the actual index sequence. Remaining columns must follow strict ontology: see 3rd tab of https://www.ncbi.nlm.nih.gov/core/assets/sra/files/SRA_metadata_acc_example.xlsx for controlled vocabulary and term definitions.",
            patterns: ["*.txt", "*.tsv"]
            category: "required"
        }
        sample_rename_map: {
            description: "If 'samples' need to be renamed, provide a two-column tsv that contains at least the following columns: internal_id, external_id. All samples will be renamed prior to analysis. Any samples described in the samplesheets that are not present in sample_rename_map will be unaltered. If this is omitted, no samples will be renamed.",
            patterns: ["*.txt", "*.tsv"]
            category: "advanced"
        }
        biosample_map: {
            description: "A two-column tsv that contains at least the following columns: sample_name, accession. sample_name refers to the external sample id, accession is the NCBI BioSample accession number (SAMNxxx). If this file is omitted, SRA submission prep will be skipped.",
            patterns: ["*.txt", "*.tsv"]
            category: "advanced"
        }
        spikein_db: {
            description: "Archeal DNA that is used to track potential contamination within sequencing plate processing."
            category: "advanced"
        }
        read_structure: { 
            description: "File that details how the bases should be organized into logical reads."
            category: "advanced"
        }
        sort_reads: {
            description: "Barcode/index information organized into barcode files to sort data to separate files for each sample."
            category: "advanced"
        }
        bmtaggerDbs: {
            description: "Tool that can discriminate between human and bacterial reads and other reads by using short fragments. Databases must be provided to onset depletion.Sequences in fasta format will be indexed on the fly, pre-bmtagger-indexed databases may be provided as tarballs."
            category: "advanced"
        }
        cleaned_bams_tiny: {
            description: "BAM files that did not meet the sufficient amount of minimum reads mapped."
            category: "other"
        }
        cleaned_bam_uris: {
            description: "BAM files' URIs tags to clearly define which assembly has been used "
            category: "other"
        }
        cleaned_reads_unaligned_bams: {
            description: "Unaligned reads without PCR duplicates in BAM format."
            category: "other"
        }
        demux_commonBarcodes: {
            description: "a TSV report of all barcode counts, in descending order."
            category: "other"
        }
        demux_metrics: { 
            description: "Output ExtractIlluminaBarcodes metrics file. Default is to dump to a temp file."
            category: "other"
        }
        multiqc_report_cleaned: {
            description: "Aggregate results from QC analyses across many samples into a single report."
            category: "other"
        }
        multiqc_report: {
            description: "Aggregation of QC metrics for all indexed samples."
            category: "other"
        }
        raw_reads_unaligned_bams : { 
            description: "Unaligned reads in BAM format."
            category: "other"
        }
        read_counts_raw: {
            description: "The number of reads that aligned to gene."
            category: "other"
        }
    }

    #### demux each lane (rename samples if requested)
    scatter(lane_sheet in zip(range(length(samplesheets)), samplesheets)) {
        call demux.samplesheet_rename_ids {
            input:
                old_sheet  = lane_sheet.right,
                rename_map = sample_rename_map
        }
        call demux.illumina_demux {
            input:
                flowcell_tgz  = flowcell_tgz,
                lane          = lane_sheet.left + 1,
                samplesheet   = samplesheet_rename_ids.new_sheet,
                readStructure = read_structure,
                sort_reads    = sort_reads
        }
        call demux.map_map_setdefault as meta_default_sample {
            input:
                map_map_json = illumina_demux.meta_by_sample_json,
                sub_keys     = default_sample_keys
        }
        call demux.map_map_setdefault as meta_default_filename {
            input:
                map_map_json = illumina_demux.meta_by_filename_json,
                sub_keys     = default_filename_keys
        }
    }
    call demux.merge_maps as meta_sample {
        input: maps_jsons = meta_default_sample.out_json
    }
    call demux.merge_maps as meta_filename {
        input: maps_jsons = meta_default_filename.out_json
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
        if (deplete.depletion_read_count_post >= min_reads_per_bam) {
            File cleaned_bam_passing = deplete.cleaned_bam
        }
        if (deplete.depletion_read_count_post < min_reads_per_bam) {
            File empty_bam = raw_reads
        }
    }

    #### SRA submission prep
    if(defined(biosample_map)) {
        call ncbi.sra_meta_prep {
            input:
                cleaned_bam_filepaths = select_all(cleaned_bam_passing),
                biosample_map         = select_first([biosample_map]),
                library_metadata      = samplesheet_rename_ids.new_sheet,
                platform              = "ILLUMINA",
                paired                = (illumina_demux.run_info[0]['indexes'] == '2'),

                out_name              = "sra_metadata-~{illumina_demux.run_info[0]['run_id']}.tsv",
                instrument_model      = select_first(flatten([[instrument_model_user_specified],[illumina_demux.run_info[0]['sequencer_model']]])),
                title                 = select_first([sra_title])
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
        Array[File] raw_reads_unaligned_bams                 = flatten(illumina_demux.raw_reads_unaligned_bams)
        Array[Int]  read_counts_raw                          = deplete.depletion_read_count_pre
        
        Map[String,Map[String,String]] meta_by_filename      = meta_filename.merged
        Map[String,Map[String,String]] meta_by_sample        = meta_sample.merged
        File                           meta_by_filename_json = meta_filename.merged_json
        File                           meta_by_sample_json   = meta_sample.merged_json
        
        Array[File] cleaned_reads_unaligned_bams             = select_all(cleaned_bam_passing)
        Array[File] cleaned_bams_tiny                        = select_all(empty_bam)
        Array[Int]  read_counts_depleted                     = deplete.depletion_read_count_post
        
        File?       sra_metadata                             = sra_meta_prep.sra_metadata
        File?       cleaned_bam_uris                         = sra_meta_prep.cleaned_bam_uris
        
        Array[File] demux_metrics                            = illumina_demux.metrics
        Array[File] demux_commonBarcodes                     = illumina_demux.commonBarcodes
        Array[File] demux_outlierBarcodes                    = illumina_demux.outlierBarcodes
        
        File        multiqc_report_raw                       = multiqc_raw.multiqc_report
        File        multiqc_report_cleaned                   = multiqc_cleaned.multiqc_report
        File        spikein_counts                           = spike_summary.count_summary
        
        String      instrument_model_inferred                = select_first(flatten([[instrument_model_user_specified],[illumina_demux.run_info[0]['sequencer_model']]]))

        String             run_date                          = illumina_demux.run_info[0]['run_start_date']
        Map[String,String] run_info                          = illumina_demux.run_info[0]
        File               run_info_json                     = illumina_demux.run_info_json[0]
        String             run_id                            = illumina_demux.run_info[0]['run_id']
        
        String      demux_viral_core_version                 = illumina_demux.viralngs_version[0]
    }
}
