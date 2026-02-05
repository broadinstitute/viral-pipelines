version 1.0

#DX_SKIP_WORKFLOW

import "../tasks/tasks_demux.wdl" as demux
import "../tasks/tasks_reports.wdl" as reports

workflow load_illumina_fastqs {
  meta {
    description: "Process DRAGEN FASTQ files into per-sample BAM files with demultiplexing support."
    author: "Broad Viral Genomics"
    email:  "viral-ngs@broadinstitute.org"
    allowNestedInputs: true
  }

  input {
    Array[File] fastq_files   # FASTQ files to demultiplex
    File        samplesheet   # TSV samplesheet with barcodes
    File        runinfo_xml   # Illumina RunInfo.xml

    Int         demux_max_cpu_splitcode    = 64   # CPU cap for 3-barcode samples (splitcode)
    Int         demux_max_cpu_no_splitcode = 16   # CPU cap for 2-barcode samples (samtools import)

    Boolean     run_fastqc = true   # Run FastQC/MultiQC reports (can be disabled for speed)
  }

  # Step 1: Group FASTQs into R1/R2 pairs (convert Files to Strings to avoid localization)
  call demux.group_fastq_pairs {
    input:
      fastq_uris = fastq_files
  }

  # Step 2: Extract run-level metadata (flowcell info, tile counts, etc.)
  call demux.get_illumina_run_metadata {
    input:
      runinfo_xml = runinfo_xml
  }

  # Step 2b: Check if samplesheet has barcode_3 (determines demux CPU allocation)
  call demux.check_for_barcode3 {
    input:
      samplesheet = samplesheet
  }

  # Step 3: Demux each FASTQ pair in parallel
  scatter (fastq_pair in group_fastq_pairs.paired_fastqs) {
    # Create a null File? variable for single-end FASTQs
    if (false) {
      File null_file = fastq_pair[0]
    }

    call demux.demux_fastqs {
      input:
        fastq_r1    = fastq_pair[0],
        fastq_r2    = if length(fastq_pair) > 1 then fastq_pair[1] else null_file,
        samplesheet = samplesheet,
        runinfo_xml = runinfo_xml,
        max_cpu     = if check_for_barcode3.has_barcode3 then demux_max_cpu_splitcode else demux_max_cpu_no_splitcode
    }
  }

  # Step 4: Merge demux metrics from all FASTQ pairs
  call demux.merge_demux_metrics {
    input:
      metrics_files = demux_fastqs.demux_metrics
  }

  # Step 4b: Merge sample metadata from all FASTQ pairs
  call demux.merge_sample_metadata {
    input:
      meta_by_sample_jsons   = demux_fastqs.meta_by_sample_json,
      meta_by_filename_jsons = demux_fastqs.meta_by_filename_json
  }

  # Flatten BAMs for convenience
  Array[File] raw_bams = flatten(demux_fastqs.output_bams)

  # Step 5: Aggregate QC reports with MultiQC (FastQC + MultiQC in one step)
  if (run_fastqc) {
    call reports.multiqc_from_bams as multiqc {
      input:
        input_bams = raw_bams
    }
  }

  output {
    # BAM outputs (flattened)
    Array[File] raw_reads_unaligned_bams = raw_bams
    Array[Int]  read_counts_raw          = flatten(demux_fastqs.read_counts)

    # QC outputs (FastQC HTML files from multiqc_from_bams)
    Array[File] fastqcs = select_first([multiqc.fastqc_html, []])

    # MultiQC aggregated report
    File? multiqc_report      = multiqc.multiqc_report
    File? multiqc_data_tar_gz = multiqc.multiqc_data_dir_tarball

    # Demux metrics
    File demux_metrics = merge_demux_metrics.merged_metrics

    # Metadata outputs (File only - Map types not supported by womtool)
    File   meta_by_sample_json   = merge_sample_metadata.merged_meta_by_sample
    File   meta_by_filename_json = merge_sample_metadata.merged_meta_by_filename
    File   run_info_json         = get_illumina_run_metadata.runinfo_json

    # Run info as Map (simple Map is supported)
    Map[String,String] run_info  = get_illumina_run_metadata.run_info

    String viralngs_version = get_illumina_run_metadata.viralngs_version
  }
}
