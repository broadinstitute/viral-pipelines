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
    Array[File]? bmtaggerDbs
    Array[File]? blastDbs
    Array[File]? bwaDbs

    Boolean      insert_demux_outputs_into_terra_tables = false
    Int          min_reads_per_bam = 100

    Array[File]  biosample_map_tsvs = []
    String?      instrument_model_user_specified
    String?      sra_title
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

  # Step 3: Extract run metadata
  call demux.get_illumina_run_metadata {
    input:
      samplesheet = samplesheet_rename_ids.new_sheet,
      runinfo_xml = runinfo_xml
  }

  # Step 4: Set default metadata keys
  call demux.map_map_setdefault as meta_default_sample {
    input:
      map_map_json = get_illumina_run_metadata.meta_by_sample_json,
      sub_keys     = default_sample_keys
  }
  call demux.map_map_setdefault as meta_default_filename {
    input:
      map_map_json = get_illumina_run_metadata.meta_by_filename_json,
      sub_keys     = default_filename_keys
  }

  # Step 5: Demux each FASTQ pair in parallel
  scatter (fastq_pair in group_fastq_pairs.paired_fastqs) {
    if (false) { File null_file = fastq_pair[0] }

    call demux.demux_fastqs {
      input:
        fastq_r1    = fastq_pair[0],
        fastq_r2    = if length(fastq_pair) > 1 then fastq_pair[1] else null_file,
        samplesheet = samplesheet_rename_ids.new_sheet,
        runinfo_xml = runinfo_xml
    }
  }

  Array[File] raw_bams = flatten(demux_fastqs.output_bams)

  # Step 6: Spike-in counting and depletion for all BAMs
  scatter (raw_reads in raw_bams) {
    # Spike-in counting (if spikein_db provided)
    if (defined(spikein_db)) {
      call reports.align_and_count as spikein {
        input:
          reads_bam = raw_reads,
          ref_db = select_first([spikein_db])
      }
    }
    Int reads_total_count = select_first([spikein.reads_total, 0])

    # Depletion (if any depletion db provided)
    if (length(flatten(select_all([bmtaggerDbs, blastDbs, bwaDbs]))) > 0) {
      call taxon_filter.deplete_taxa as deplete {
        input:
          raw_reads_unmapped_bam = raw_reads,
          bmtaggerDbs = bmtaggerDbs,
          blastDbs = blastDbs,
          bwaDbs = bwaDbs
      }
    }

    Int  read_count_post_depletion = select_first([deplete.depletion_read_count_post, reads_total_count])
    File cleaned_bam = select_first([deplete.cleaned_bam, raw_reads])

    if (read_count_post_depletion >= min_reads_per_bam) {
      File cleaned_bam_passing = cleaned_bam
    }
    if (read_count_post_depletion < min_reads_per_bam) {
      File empty_bam = raw_reads
    }

    Pair[String,Int] count_raw = (basename(raw_reads, '.bam'), reads_total_count)
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

  # Step 9: MultiQC reports
  call reports.MultiQC as multiqc_raw {
    input:
      input_files = flatten(demux_fastqs.fastqc_zip),
      file_name   = "multiqc-raw.html"
  }

  if (length(flatten(select_all([bmtaggerDbs, blastDbs, bwaDbs]))) > 0) {
    call reports.MultiQC as multiqc_cleaned {
      input:
        input_files = select_all(deplete.cleaned_fastqc_zip),
        file_name   = "multiqc-cleaned.html"
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
    Array[Int]  read_counts_raw          = reads_total_count

    # Metadata outputs (with defaults applied)
    Map[String,Map[String,String]] meta_by_sample   = read_json(meta_default_sample.out_json)
    Map[String,Map[String,String]] meta_by_filename = read_json(meta_default_filename.out_json)
    File                           meta_by_sample_json   = meta_default_sample.out_json
    File                           meta_by_filename_json = meta_default_filename.out_json

    # Cleaned BAM outputs
    Array[File] cleaned_reads_unaligned_bams = select_all(cleaned_bam_passing)
    Array[File] cleaned_bams_tiny            = select_all(empty_bam)
    Array[Int]  read_counts_depleted         = read_count_post_depletion

    # SRA outputs
    File? sra_metadata     = sra_meta_prep.sra_metadata
    File? cleaned_bam_uris = sra_meta_prep.cleaned_bam_uris

    # QC outputs
    Array[File] fastqc_html = flatten(demux_fastqs.fastqc_html)
    Array[File] fastqc_zip  = flatten(demux_fastqs.fastqc_zip)

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
