version 1.0

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
  }

  # Step 1: Group FASTQs into R1/R2 pairs (convert Files to Strings to avoid localization)
  call demux.group_fastq_pairs {
    input:
      fastq_uris = fastq_files
  }

  # Step 2: Extract run metadata (once for entire run)
  call demux.get_illumina_run_metadata {
    input:
      samplesheet = samplesheet,
      runinfo_xml = runinfo_xml
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
        runinfo_xml = runinfo_xml
    }
  }

  # Step 4: Aggregate QC reports with MultiQC
  call reports.MultiQC {
    input:
      input_files = flatten(demux_fastqs.fastqc_zip)
  }

  output {
    # BAM outputs (flattened)
    Array[File] raw_reads_unaligned_bams = flatten(demux_fastqs.output_bams)
    Array[Int]  read_counts_raw          = flatten(demux_fastqs.read_counts)

    # QC outputs (flattened)
    Array[File] fastqc_html = flatten(demux_fastqs.fastqc_html)
    Array[File] fastqc_zip  = flatten(demux_fastqs.fastqc_zip)

    # MultiQC aggregated report
    File multiqc_report      = MultiQC.multiqc_report
    File multiqc_data_tar_gz = MultiQC.multiqc_data_dir_tarball

    # Metadata outputs
    Map[String,Map[String,String]] meta_by_sample   = get_illumina_run_metadata.meta_by_sample
    Map[String,Map[String,String]] meta_by_filename = get_illumina_run_metadata.meta_by_filename
    Map[String,String]             run_info         = get_illumina_run_metadata.run_info

    File   meta_by_sample_json   = get_illumina_run_metadata.meta_by_sample_json
    File   meta_by_filename_json = get_illumina_run_metadata.meta_by_filename_json
    File   run_info_json         = get_illumina_run_metadata.runinfo_json

    String viralngs_version = get_illumina_run_metadata.viralngs_version
  }
}
