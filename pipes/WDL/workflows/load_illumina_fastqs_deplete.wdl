version 1.0

#DX_SKIP_WORKFLOW

import "../tasks/tasks_demux.wdl" as demux
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_terra.wdl" as terra
import "../tasks/tasks_utils.wdl" as utils

workflow load_illumina_fastqs_deplete {
  meta {
    description: "Process DRAGEN FASTQ files into per-sample BAM files with demultiplexing, depletion, and SRA submission prep."
    author: "Broad Viral Genomics"
    email:  "viral-ngs@broadinstitute.org"
    allowNestedInputs: true
  }

  input {
    Array[File] fastq_files
    File        samplesheet
    File        runinfo_xml

    File?        sample_rename_map

    Array[String] default_sample_keys = ["amplicon_set", "control", "batch_lib", "viral_ct"]
    Array[String] default_filename_keys = ["spike_in", "batch_lib"]

    File?        spikein_db
    Array[File]? minimapDbs
    Array[File]? bmtaggerDbs
    Array[File]? blastDbs
    Array[File]? bwaDbs

    Boolean      insert_demux_outputs_into_terra_tables = false
    Int          min_reads_per_bam = 100

    Array[File]  biosample_map_tsvs = []
    String?      instrument_model_user_specified
    String?      sra_title

    Int          demux_max_cpu_splitcode    = 64   # CPU cap for 3-barcode samples (splitcode)
    Int          demux_max_cpu_no_splitcode = 16   # CPU cap for 2-barcode samples (samtools import)
  }

  # Step 1: Rename samples in samplesheet (if sample_rename_map provided)
  call demux.samplesheet_rename_ids {
    input:
      old_sheet  = samplesheet,
      rename_map = sample_rename_map
  }

  # Step 2: Group FASTQs into R1/R2 pairs
  call demux.group_fastq_pairs {
    input:
      fastq_uris = fastq_files
  }

  # Step 3: Extract run-level metadata (flowcell info, tile counts, etc.)
  call demux.get_illumina_run_metadata {
    input:
      runinfo_xml = runinfo_xml
  }

  # Step 3b: Check if samplesheet has barcode_3 (determines demux CPU allocation)
  call demux.check_for_barcode3 {
    input:
      samplesheet = samplesheet_rename_ids.new_sheet
  }

  # Step 5: Demux each FASTQ pair in parallel
  scatter (fastq_pair in group_fastq_pairs.paired_fastqs) {
    if (false) { File null_file = fastq_pair[0] }

    call demux.demux_fastqs {
      input:
        fastq_r1    = fastq_pair[0],
        fastq_r2    = if length(fastq_pair) > 1 then fastq_pair[1] else null_file,
        samplesheet = samplesheet_rename_ids.new_sheet,
        runinfo_xml = runinfo_xml,
        max_cpu     = if check_for_barcode3.has_barcode3 then demux_max_cpu_splitcode else demux_max_cpu_no_splitcode
    }
  }

  Array[File] raw_bams = flatten(demux_fastqs.output_bams)
  Array[Int]  raw_read_counts = flatten(demux_fastqs.read_counts)

  # Step 5.5: Merge demux metrics from all FASTQ pairs
  call demux.merge_demux_metrics {
    input:
      metrics_files = demux_fastqs.demux_metrics
  }

  # Step 5.6: Merge sample metadata from all FASTQ pairs
  call demux.merge_sample_metadata {
    input:
      meta_by_sample_jsons   = demux_fastqs.meta_by_sample_json,
      meta_by_filename_jsons = demux_fastqs.meta_by_filename_json
  }

  # Step 5.7: Set default metadata keys
  call demux.map_map_setdefault as meta_default_sample {
    input:
      map_map_json = merge_sample_metadata.merged_meta_by_sample,
      sub_keys     = default_sample_keys
  }
  call demux.map_map_setdefault as meta_default_filename {
    input:
      map_map_json = merge_sample_metadata.merged_meta_by_filename,
      sub_keys     = default_filename_keys
  }

  # Step 6: Spike-in counting and depletion for all BAMs
  scatter (idx in range(length(raw_bams))) {
    File raw_reads = raw_bams[idx]
    Int  raw_count = raw_read_counts[idx]

    # Spike-in counting (if spikein_db provided)
    if (defined(spikein_db)) {
      call reports.align_and_count as spikein {
        input:
          reads_bam = raw_reads,
          ref_db = select_first([spikein_db])
      }
    }

    # Depletion (if any depletion db provided)
    if (length(flatten(select_all([minimapDbs, bmtaggerDbs, blastDbs, bwaDbs]))) > 0) {
      call taxon_filter.deplete_taxa as deplete {
        input:
          raw_reads_unmapped_bam = raw_reads,
          minimapDbs = minimapDbs,
          bmtaggerDbs = bmtaggerDbs,
          blastDbs = blastDbs,
          bwaDbs = bwaDbs
      }
    }

    Int  read_count_post_depletion = select_first([deplete.depletion_read_count_post, raw_count])
    File cleaned_bam = select_first([deplete.cleaned_bam, raw_reads])

    if (read_count_post_depletion >= min_reads_per_bam) {
      File cleaned_bam_passing = cleaned_bam
    }
    if (read_count_post_depletion < min_reads_per_bam) {
      File empty_bam = raw_reads
    }

    Pair[String,Int] count_raw = (basename(raw_reads, '.bam'), raw_count)
    Pair[String,Int] count_cleaned = (basename(raw_reads, '.bam'), read_count_post_depletion)
  }

  # Step 7: Terra table insertion (if enabled and running on Terra)
  if (insert_demux_outputs_into_terra_tables) {
    call terra.check_terra_env

    if (check_terra_env.is_running_on_terra) {
      call terra.create_or_update_sample_tables {
        input:
          flowcell_run_id          = get_illumina_run_metadata.run_info['run_id'],
          workspace_name           = check_terra_env.workspace_name,
          workspace_namespace      = check_terra_env.workspace_namespace,
          raw_reads_unaligned_bams     = raw_bams,
          cleaned_reads_unaligned_bams = select_all(cleaned_bam_passing),
          meta_by_filename_json        = meta_default_filename.out_json,
          read_counts_raw_json         = write_json(count_raw),
          read_counts_cleaned_json     = write_json(count_cleaned)
      }
    }
  }

  # Step 8: BioSample/SRA prep (if biosample_map_tsvs provided)
  if (length(biosample_map_tsvs) > 0) {
    # Merge biosample tsvs if more than one
    if (length(biosample_map_tsvs) > 1) {
      call utils.tsv_join as biosample_map_tsv_join {
        input:
          input_tsvs   = biosample_map_tsvs,
          id_col       = 'accession',
          out_suffix   = ".tsv",
          out_basename = "biosample-attributes-merged"
      }
    }
    File biosample_map_tsv = select_first(flatten([[biosample_map_tsv_join.out_tsv], biosample_map_tsvs]))

    call ncbi.biosample_to_table {
      input:
        biosample_attributes_tsv = biosample_map_tsv,
        raw_bam_filepaths        = raw_bams,
        demux_meta_json          = meta_default_filename.out_json
    }

    call ncbi.sra_meta_prep {
      input:
        cleaned_bam_filepaths = select_all(cleaned_bam_passing),
        biosample_map         = biosample_map_tsv,
        library_metadata      = [samplesheet_rename_ids.new_sheet],
        platform              = "ILLUMINA",
        paired                = (get_illumina_run_metadata.run_info['indexes'] == '2'),
        out_name              = "sra_metadata-~{get_illumina_run_metadata.run_info['run_id']}.tsv",
        instrument_model      = select_first(flatten([[instrument_model_user_specified], [get_illumina_run_metadata.run_info['sequencer_model']]])),
        title                 = select_first([sra_title, "Viral sequencing"])
    }

    if (insert_demux_outputs_into_terra_tables && select_first([check_terra_env.is_running_on_terra, false])) {
      call terra.upload_entities_tsv as terra_load_biosample_data {
        input:
          workspace_name = select_first([check_terra_env.workspace_name]),
          terra_project  = select_first([check_terra_env.workspace_namespace]),
          tsv_file       = biosample_to_table.sample_meta_tsv
      }
    }
  }

  # Step 9: MultiQC reports (FastQC + MultiQC in one step)
  call reports.multiqc_from_bams as multiqc_raw {
    input:
      input_bams   = raw_bams,
      out_basename = "multiqc-raw"
  }

  if (length(flatten(select_all([minimapDbs, bmtaggerDbs, blastDbs, bwaDbs]))) > 0) {
    call reports.multiqc_from_bams as multiqc_cleaned {
      input:
        input_bams   = select_all(cleaned_bam_passing),
        out_basename = "multiqc-cleaned"
    }
  }

  # Step 10: Spike-in summary (if spikein_db provided)
  if (defined(spikein_db)) {
    call reports.align_and_count_summary as spike_summary {
      input:
        counts_txt = select_all(spikein.report)
    }
  }

  output {
    # Raw BAM outputs
    Array[File] raw_reads_unaligned_bams = raw_bams
    Array[Int]  read_counts_raw          = raw_read_counts

    # Metadata outputs (File only - Map types not supported by womtool)
    File meta_by_sample_json   = meta_default_sample.out_json
    File meta_by_filename_json = meta_default_filename.out_json

    # Cleaned BAM outputs
    Array[File] cleaned_reads_unaligned_bams = select_all(cleaned_bam_passing)
    Array[File] cleaned_bams_tiny            = select_all(empty_bam)
    Array[Int]  read_counts_depleted         = read_count_post_depletion

    # SRA outputs
    File? sra_metadata     = sra_meta_prep.sra_metadata
    File? cleaned_bam_uris = sra_meta_prep.cleaned_bam_uris

    # QC outputs (FastQC HTML files from multiqc_from_bams)
    Array[File] fastqcs_raw     = multiqc_raw.fastqc_html
    Array[File] fastqcs_cleaned = select_first([multiqc_cleaned.fastqc_html, []])

    # Demux metrics
    File demux_metrics = merge_demux_metrics.merged_metrics

    # MultiQC reports
    File  multiqc_report_raw     = multiqc_raw.multiqc_report
    File  multiqc_report_cleaned = select_first([multiqc_cleaned.multiqc_report, multiqc_raw.multiqc_report])

    # Spike-in outputs
    File? spikein_counts = spike_summary.count_summary

    # Run info
    Map[String,String] run_info      = get_illumina_run_metadata.run_info
    File               run_info_json = get_illumina_run_metadata.runinfo_json
    String             viralngs_version = get_illumina_run_metadata.viralngs_version

    # Terra outputs
    File? terra_library_table      = create_or_update_sample_tables.library_metadata_tsv
    File? terra_sample_library_map = create_or_update_sample_tables.sample_membership_tsv
    File? terra_sample_metadata    = biosample_to_table.sample_meta_tsv

    # Instrument model
    String instrument_model_inferred = select_first(flatten([[instrument_model_user_specified], [get_illumina_run_metadata.run_info['sequencer_model']]]))
  }
}
