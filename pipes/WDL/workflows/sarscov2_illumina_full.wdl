version 1.0

import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_sarscov2.wdl" as sarscov2
import "../tasks/tasks_terra.wdl" as terra
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_utils.wdl" as utils

import "demux_deplete.wdl"
import "assemble_refbased.wdl"
import "sarscov2_lineages.wdl"

workflow sarscov2_illumina_full {
    meta {
        description: "Full SARS-CoV-2 analysis workflow starting from raw Illumina flowcell (tar.gz) and metadata and performing assembly, spike-in analysis, qc, lineage assignment, and packaging for data release."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    parameter_meta {
        reference_fasta: {
            description: "Reference genome to align reads to.",
            patterns: ["*.fasta"]
        }
        amplicon_bed_prefix: {
            description: "amplicon primers to trim in reference coordinate space (0-based BED format)",
            patterns: ["*.bed"]
        }

        biosample_attributes: {
          description: "A post-submission attributes file from NCBI BioSample, which is available at https://submit.ncbi.nlm.nih.gov/subs/ and clicking on 'Download attributes file with BioSample accessions'. The 'sample_name' column must match the external_ids used in sample_rename_map (or internal ids if sample_rename_map is omitted).",
          patterns: ["*.txt", "*.tsv"]
        }
        author_list: {
          description: "A string containing a space-delimited list with of author surnames separated by first name and (optional) middle initial. Ex. 'Lastname,Firstname, Last-hypenated,First,M., Last,F.'"
        }
        author_sbt_defaults_yaml: {
          description: "A YAML file with default values to use for the submitter, submitter affiliation, and author affiliation. Optionally including authors at the start and end of the author_list. Example: gs://pathogen-public-dbs/other-related/default_sbt_values.yaml",
          patterns: ["*.yaml","*.yml"]
        }
        author_sbt_j2_template: {
          description: "A jinja2-format template for the sbt file expected by NCBI. Example: gs://pathogen-public-dbs/other-related/author_template.sbt.j2"
        }
    }

    input {
        File          flowcell_tgz
        File          reference_fasta
        String        amplicon_bed_prefix

        String?       author_list # of the form "Lastname,A.B., Lastname,C.,"; optional alternative to names in author_sbt_defaults_yaml
        File          author_sbt_defaults_yaml # defaults to fill in for author_sbt file (including both author and non-author fields)
        File          author_sbt_j2_template

        Array[File]   biosample_attributes
        String?       instrument_model
        String        sra_title

        Int           min_genome_bases = 24000
        Int           max_vadr_alerts = 0


        String?       workspace_name
        String?       terra_project
        File?         collab_ids_tsv

        String?       gcs_out_metrics
        String?       gcs_out_reporting
        String?       gcs_out_cdc
        String?       gcs_out_sra
    }
    Int     taxid         = 2697049
    String  gisaid_prefix = 'hCoV-19/'

    # merge biosample attributes tables
    call utils.tsv_join as biosample_merge {
        input:
            input_tsvs   = biosample_attributes,
            id_col       = 'accession',
            out_basename = "biosample_attributes-merged"
    }
    call utils.fetch_col_from_tsv as accessioned_samples {
      input:
        tsv = biosample_merge.out_tsv,
        col = 'sample_name'
    }

    ### demux, deplete, SRA submission prep, fastqc/multiqc
    call demux_deplete.demux_deplete {
        input:
            flowcell_tgz                    = flowcell_tgz,
            biosample_map                   = biosample_merge.out_tsv,
            instrument_model_user_specified = instrument_model,
            sra_title                       = sra_title
    }
    String  flowcell_id = demux_deplete.run_id

    ### gather data by biosample
    call read_utils.group_bams_by_sample {
        input:
            bam_filepaths = demux_deplete.cleaned_reads_unaligned_bams
    }

    ### assembly and analyses per biosample
    scatter(name_reads in zip(group_bams_by_sample.sample_names, group_bams_by_sample.grouped_bam_filepaths)) {
        Boolean ampseq   = (demux_deplete.meta_by_sample[name_reads.left]["amplicon_set"] != "")
        String orig_name = demux_deplete.meta_by_sample[name_reads.left]["sample_original"]

        # assemble genome
        if (ampseq) {
            String trim_coords_bed = amplicon_bed_prefix + demux_deplete.meta_by_sample[name_reads.left]["amplicon_set"] + ".bed"
        }
        call assemble_refbased.assemble_refbased {
            input:
                reads_unmapped_bams = name_reads.right,
                reference_fasta     = reference_fasta,
                sample_name         = name_reads.left,
                aligner             = "minimap2",
                skip_mark_dupes     = ampseq,
                trim_coords_bed     = trim_coords_bed,
                min_coverage        = if ampseq then 20 else 3
        }

        # log controls
        if (demux_deplete.meta_by_sample[name_reads.left]["control"] == 'NTC') {
            Int ntc_bases = assemble_refbased.assembly_length_unambiguous
        }

        # grab biosample metadata
        call utils.fetch_row_from_tsv as biosample {
            input:
                tsv              = biosample_merge.out_tsv,
                idx_col          = "sample_name",
                idx_val          = orig_name,
                set_default_keys = ["collection_date", "bioproject_accession", "accession", "collected_by", "geo_loc_name", "host_subject_id", "host_age", "host_sex", "purpose_of_sequencing"]
        }

        # for genomes that somewhat assemble
        if (assemble_refbased.assembly_length_unambiguous >= min_genome_bases) {
            call ncbi.rename_fasta_header {
              input:
                genome_fasta = assemble_refbased.assembly_fasta,
                new_name     = orig_name
            }

            File passing_assemblies     = rename_fasta_header.renamed_fasta
            String passing_assembly_ids = orig_name
            Array[String] assembly_cmt  = [orig_name, "Broad viral-ngs v. " + demux_deplete.demux_viral_core_version, assemble_refbased.assembly_mean_coverage, demux_deplete.instrument_model_inferred]

            # lineage assignment
            call sarscov2_lineages.sarscov2_lineages {
                input:
                    genome_fasta = passing_assemblies
            }

            # VADR annotation & QC
            call ncbi.vadr {
              input:
                genome_fasta = passing_assemblies
            }
            if (vadr.num_alerts<=max_vadr_alerts) {
              File submittable_genomes = passing_assemblies
              String submittable_id    = orig_name
            }
            if (vadr.num_alerts>max_vadr_alerts) {
              String failed_annotation_id = orig_name
            }
        }
        if (assemble_refbased.assembly_length_unambiguous < min_genome_bases) {
            String failed_assembly_id = orig_name
        }

        Array[String] assembly_tsv_row = [
            orig_name,
            name_reads.left,
            biosample.map["accession"],
            flowcell_id,
            demux_deplete.run_date,
            biosample.map["collection_date"],
            biosample.map["geo_loc_name"],
            biosample.map["host_subject_id"],
            assemble_refbased.assembly_length_unambiguous,
            assemble_refbased.assembly_mean_coverage,
            select_first([sarscov2_lineages.pango_lineage, ""]),
            select_first([sarscov2_lineages.nextclade_clade, ""]),
            select_first([sarscov2_lineages.nextclade_aa_subs, ""]),
            select_first([sarscov2_lineages.nextclade_aa_dels, ""]),
            select_first([sarscov2_lineages.pangolin_docker, ""]),
            assemble_refbased.dist_to_ref_snps,
            assemble_refbased.dist_to_ref_indels,
            select_first([vadr.num_alerts, ""]),
            assemble_refbased.assembly_fasta,
            assemble_refbased.align_to_ref_merged_coverage_plot,
            assemble_refbased.align_to_ref_merged_aligned_trimmed_only_bam,
            assemble_refbased.replicate_discordant_vcf,
            assemble_refbased.align_to_ref_variants_vcf_gz,
            select_first([sarscov2_lineages.nextclade_tsv, ""]),
            select_first([sarscov2_lineages.nextclade_json, ""]),
            select_first([sarscov2_lineages.pangolin_csv, ""]),
            select_first([vadr.outputs_tgz, ""]),
            demux_deplete.meta_by_sample[name_reads.left]["amplicon_set"],
            assemble_refbased.replicate_concordant_sites,
            assemble_refbased.replicate_discordant_snps,
            assemble_refbased.replicate_discordant_indels,
            assemble_refbased.num_read_groups,
            assemble_refbased.num_libraries,
            assemble_refbased.align_to_ref_merged_reads_aligned,
            assemble_refbased.align_to_ref_merged_bases_aligned,
            select_first([vadr.alerts_list, ""]),
            biosample.map["purpose_of_sequencing"],
            biosample.map["collected_by"],
            biosample.map["bioproject_accession"],
            biosample.map["host_age"],
            biosample.map["host_sex"],
            "",
            demux_deplete.meta_by_sample[name_reads.left]["viral_ct"]
        ]
    }
    Array[String] assembly_tsv_header = [
        'sample', 'sample_sanitized', 'biosample_accession', 'flowcell_id', 'run_date', 'collection_date', 'geo_loc_name', 'host_subject_id',
        'assembly_length_unambiguous', 'assembly_mean_coverage',
        'pango_lineage', 'nextclade_clade', 'nextclade_aa_subs', 'nextclade_aa_dels', 'pangolin_version',
        'dist_to_ref_snps', 'dist_to_ref_indels', 'vadr_num_alerts',
        'assembly_fasta', 'coverage_plot', 'aligned_bam',
        'replicate_discordant_vcf', 'variants_from_ref_vcf',
        'nextclade_tsv', 'nextclade_json', 'pangolin_csv', 'vadr_tgz',
        'amplicon_set',
        'replicate_concordant_sites', 'replicate_discordant_snps', 'replicate_discordant_indels', 'num_read_groups', 'num_libraries',
        'align_to_ref_merged_reads_aligned', 'align_to_ref_merged_bases_aligned',
        'vadr_alerts', 'purpose_of_sequencing', 'collected_by', 'bioproject_accession',
        'age', 'sex', 'zip',
        'Ct'
        ]

    ### summary stats
    call utils.concatenate as assembly_meta_tsv {
      input:
        infiles     = [write_tsv([assembly_tsv_header]), write_tsv(assembly_tsv_row)],
        output_name = "assembly_metadata-~{flowcell_id}.tsv"
    }

    ### mark up the bad batches or lanes where NTCs assemble
    call assembly.filter_bad_ntc_batches {
      input:
        seqid_list = write_lines(select_all(passing_assembly_ids)),
        demux_meta_by_sample_json = demux_deplete.meta_by_sample_json,
        assembly_meta_tsv = assembly_meta_tsv.combined,
        ntc_min_unambig = 15000
    }

    ### QC metrics
    call read_utils.max as ntc_max {
      input:
        list = select_all(ntc_bases)
    }
    call assembly.ivar_trim_stats {
      input:
        ivar_trim_stats_tsv = write_tsv(flatten(assemble_refbased.ivar_trim_stats_tsv)),
        flowcell            = flowcell_id,
        out_basename        = "ivar_trim_stats-~{flowcell_id}"
    }
    call utils.tsv_join as picard_wgs_merge {
      input:
        input_tsvs   = assemble_refbased.picard_metrics_wgs,
        id_col       = 'sample_sanitized',
        out_basename = "picard_metrics_wgs-~{flowcell_id}"
    }
    call utils.tsv_join as picard_alignment_merge {
      input:
        input_tsvs   = assemble_refbased.picard_metrics_alignment,
        id_col       = 'sample_sanitized',
        out_basename = "picard_metrics_alignment-~{flowcell_id}"
    }
    call utils.tsv_join as picard_insertsize_merge {
      input:
        input_tsvs   = assemble_refbased.picard_metrics_insert_size,
        id_col       = 'sample_sanitized',
        out_basename = "picard_metrics_insertsize-~{flowcell_id}"
    }

    ### filter and concatenate final sets for delivery ("passing" and "submittable")
    call sarscov2.sc2_meta_final {
      # this decorates assembly_meta_tsv with collab/internal IDs, genome_status, and many other columns
      input:
        assembly_stats_tsv = assembly_meta_tsv.combined,
        collab_ids_tsv = collab_ids_tsv,
        drop_file_cols = true,
        min_unambig = min_genome_bases,
        genome_status_json = filter_bad_ntc_batches.fail_meta_json
    }
    call utils.concatenate as passing_cat_prefilter {
      # this emits a fasta of only genomes that pass min_unambig
      input:
        infiles     = select_all(passing_assemblies),
        output_name = "assemblies_passing-~{flowcell_id}.prefilter.fasta"
    }
    call nextstrain.filter_sequences_to_list as passing_ntc {
      # this drops all genomes that are failed_NTC
      input:
        sequences = passing_cat_prefilter.combined,
        keep_list = [filter_bad_ntc_batches.seqids_kept]
    }
    call nextstrain.filter_sequences_to_list as passing_cat {
      # this drops all genomes that don't have BioSample accessions (e.g. control libraries)
      input:
        sequences = passing_ntc.filtered_fasta,
        keep_list = [accessioned_samples.out_txt],
        out_fname = "assemblies_passing-~{flowcell_id}.fasta"
    }
    call nextstrain.filter_sequences_to_list as submittable_filter {
      # this drops all failed_annotation (aka VADR fails)
      input:
        sequences = passing_cat.filtered_fasta,
        keep_list = [write_lines(select_all(submittable_id))]
    }

    ### prep genbank submission
    call ncbi.biosample_to_genbank {
      # this takes a BioSample attributes file and emits a Genbank Source Modifier Table
      input:
        biosample_attributes = biosample_merge.out_tsv,
        num_segments         = 1,
        taxid                = taxid,
        filter_to_ids        = submittable_filter.ids_kept
    }
    call ncbi.structured_comments {
      input:
        assembly_stats_tsv = write_tsv(flatten([[['SeqID','Assembly Method','Coverage','Sequencing Technology']],select_all(assembly_cmt)])),
        filter_to_ids      = biosample_to_genbank.sample_ids
    }
    call nextstrain.filter_sequences_to_list as submit_genomes {
      input:
        sequences = submittable_filter.filtered_fasta,
        keep_list = [biosample_to_genbank.sample_ids]
    }
    call ncbi.generate_author_sbt_file as generate_author_sbt {
        input:
            author_list   = author_list,
            defaults_yaml = author_sbt_defaults_yaml,
            j2_template   = author_sbt_j2_template
    }
    call ncbi.package_genbank_ftp_submission {
      input:
        sequences_fasta          = submit_genomes.filtered_fasta,
        source_modifier_table    = biosample_to_genbank.genbank_source_modifier_table,
        author_template_sbt      = generate_author_sbt.sbt_file,
        structured_comment_table = structured_comments.structured_comment_table,
        submission_name          = flowcell_id,
        submission_uid           = flowcell_id
    }

    ### prep gisaid submission
    call ncbi.prefix_fasta_header as prefix_gisaid {
      input:
        genome_fasta = submit_genomes.filtered_fasta,
        prefix       = gisaid_prefix,
        out_basename = "gisaid-sequences-~{flowcell_id}"
    }
    call ncbi.gisaid_meta_prep {
      input:
        source_modifier_table = biosample_to_genbank.genbank_source_modifier_table,
        structured_comments   = structured_comments.structured_comment_table,
        fasta_filename        = "gisaid-sequences-~{flowcell_id}.fasta",
        out_name              = "gisaid-meta-~{flowcell_id}.csv"
    }

    # prep nextmeta-style metadata for private nextstrain build
    call nextstrain.nextmeta_prep {
      input:
        gisaid_meta   = gisaid_meta_prep.meta_csv,
        assembly_meta = assembly_meta_tsv.combined,
        out_name      = "nextmeta-~{flowcell_id}.tsv",
        filter_to_ids = filter_bad_ntc_batches.seqids_kept
    }

    # create data tables with assembly_meta_tsv if workspace name and project provided
    if (defined(workspace_name) && defined(terra_project)) {
      call terra.upload_entities_tsv as data_tables {
        input:
          workspace_name                      = select_first([workspace_name]),
          terra_project                       = select_first([terra_project]),
          tsv_file                            = assembly_meta_tsv.combined,
          cleaned_reads_unaligned_bams_string = demux_deplete.cleaned_reads_unaligned_bams,
          meta_by_filename_json               = demux_deplete.meta_by_filename_json
      }

      call terra.download_entities_tsv {
        input:
          workspace_name   = select_first([workspace_name]),
          terra_project    = select_first([terra_project]),
          table_name       = 'assemblies',
          nop_input_string = data_tables.tables[0]
      }

      call sarscov2.sequencing_report {
        input:
            assembly_stats_tsv = download_entities_tsv.tsv_file,
            collab_ids_tsv     = collab_ids_tsv,
            max_date           = demux_deplete.run_date,
            min_unambig        = min_genome_bases
      }
    }

    # create full nextclade trees on full data set
    call sarscov2.nextclade_many_samples {
        input:
            genome_fastas = assemble_refbased.assembly_fasta,
            basename      = "nextclade-~{flowcell_id}"
    }

    # bucket deliveries
    if(defined(gcs_out_metrics)) {
        call terra.gcs_copy as gcs_metrics_dump {
            input:
              infiles        = flatten([[assembly_meta_tsv.combined, sc2_meta_final.meta_tsv, ivar_trim_stats.trim_stats_tsv, demux_deplete.multiqc_report_raw, demux_deplete.multiqc_report_cleaned, demux_deplete.spikein_counts, picard_wgs_merge.out_tsv, picard_alignment_merge.out_tsv, picard_insertsize_merge.out_tsv, nextclade_many_samples.nextclade_json, nextclade_many_samples.nextclade_tsv], demux_deplete.demux_metrics]),
              gcs_uri_prefix = "~{gcs_out_metrics}/~{flowcell_id}/"
        }
    }
    if(defined(gcs_out_reporting)) {
        call utils.tsv_to_csv as meta_final_csv {
          input:
            tsv = sc2_meta_final.meta_tsv
        }
        call terra.gcs_copy as gcs_reporting_dump {
            input:
              infiles        = [meta_final_csv.csv],
              gcs_uri_prefix = "~{gcs_out_reporting}/~{flowcell_id}/"
        }
    }
    if(defined(gcs_out_cdc)) {
        call terra.gcs_copy as gcs_cdc_dump {
            input:
                infiles        = [sc2_meta_final.meta_tsv, passing_cat.filtered_fasta],
                gcs_uri_prefix = "~{gcs_out_cdc}/~{demux_deplete.run_date}/~{flowcell_id}/"
        }
        call terra.gcs_copy as gcs_cdc_dump_reads {
            input:
                infiles        = assemble_refbased.align_to_ref_merged_aligned_trimmed_only_bam,
                gcs_uri_prefix = "~{gcs_out_cdc}/~{demux_deplete.run_date}/~{flowcell_id}/rawfiles/"
        }
    }
    if(defined(gcs_out_sra)) {
        call terra.gcs_copy as gcs_sra_dump_reads {
            input:
                infiles        = demux_deplete.cleaned_reads_unaligned_bams,
                gcs_uri_prefix = "~{gcs_out_sra}/~{flowcell_id}/"
        }
        call terra.gcs_copy as gcs_sra_dump {
            input:
                infiles        = [select_first([demux_deplete.sra_metadata])],
                gcs_uri_prefix = "~{gcs_out_sra}/"
        }
    }

    output {
        Array[File]   raw_reads_unaligned_bams      = demux_deplete.raw_reads_unaligned_bams
        Array[File]   cleaned_reads_unaligned_bams  = demux_deplete.cleaned_reads_unaligned_bams
        Array[File]   cleaned_bams_tiny             = demux_deplete.cleaned_bams_tiny
        
        File          meta_by_filename_json         = demux_deplete.meta_by_filename_json
        
        Array[Int]    read_counts_raw               = demux_deplete.read_counts_raw
        Array[Int]    read_counts_depleted          = demux_deplete.read_counts_depleted
        
        File          sra_metadata                 = select_first([demux_deplete.sra_metadata])
        File          cleaned_bam_uris             = select_first([demux_deplete.cleaned_bam_uris])
        
        Array[File]   assemblies_fasta             = assemble_refbased.assembly_fasta
        
        Int           max_ntc_bases                = ntc_max.out
        Array[String] ntc_rejected_batches         = filter_bad_ntc_batches.reject_batches
        Array[String] ntc_rejected_lanes           = filter_bad_ntc_batches.reject_lanes
        
        Array[File]   demux_metrics                = demux_deplete.demux_metrics
        Array[File]   demux_commonBarcodes         = demux_deplete.demux_commonBarcodes
        Array[File]   demux_outlierBarcodes        = demux_deplete.demux_outlierBarcodes
        
        Array[Int]    primer_trimmed_read_count    = flatten(assemble_refbased.primer_trimmed_read_count)
        Array[Float]  primer_trimmed_read_percent  = flatten(assemble_refbased.primer_trimmed_read_percent)
        File          ivar_trim_stats_html         = ivar_trim_stats.trim_stats_html
        File          ivar_trim_stats_png          = ivar_trim_stats.trim_stats_png
        File          ivar_trim_stats_tsv          = ivar_trim_stats.trim_stats_tsv
        
        File          multiqc_report_raw           = demux_deplete.multiqc_report_raw
        File          multiqc_report_cleaned       = demux_deplete.multiqc_report_cleaned
        File          spikein_counts               = demux_deplete.spikein_counts
        File          picard_metrics_wgs           = picard_wgs_merge.out_tsv
        File          picard_metrics_alignment     = picard_alignment_merge.out_tsv
        
        File          assembly_stats_tsv           = assembly_meta_tsv.combined
        File          assembly_stats_final_tsv     = sc2_meta_final.meta_tsv
        
        File          submission_zip               = package_genbank_ftp_submission.submission_zip
        File          submission_xml               = package_genbank_ftp_submission.submission_xml
        File          submit_ready                 = package_genbank_ftp_submission.submit_ready
        Array[File]   vadr_outputs                 = select_all(vadr.outputs_tgz)
        File          genbank_source_table         = biosample_to_genbank.genbank_source_modifier_table
        
        File          gisaid_fasta                 = prefix_gisaid.renamed_fasta
        File          gisaid_meta_csv              = gisaid_meta_prep.meta_csv
        
        File          genbank_fasta                = submit_genomes.filtered_fasta
        File          nextmeta_tsv                 = nextmeta_prep.nextmeta_tsv
        
        File          nextclade_all_json           = nextclade_many_samples.nextclade_json
        File          nextclade_auspice_json       = nextclade_many_samples.auspice_json
        
        File          passing_fasta                = passing_cat.filtered_fasta
        
        Array[String] assembled_ids                = select_all(passing_assembly_ids)
        Array[String] submittable_ids              = read_lines(filter_bad_ntc_batches.seqids_kept)
        Array[String] failed_assembly_ids          = select_all(failed_assembly_id)
        Array[String] failed_annotation_ids        = select_all(failed_annotation_id)
        Int           num_read_files               = length(demux_deplete.cleaned_reads_unaligned_bams)
        Int           num_assembled                = length(select_all(passing_assemblies))
        Int           num_failed_assembly          = length(select_all(failed_assembly_id))
        Int           num_submittable              = filter_bad_ntc_batches.num_kept
        Int           num_failed_annotation        = length(select_all(failed_annotation_id))
        Int           num_samples                  = length(group_bams_by_sample.sample_names)
        
        String        run_date                     = demux_deplete.run_date
        
        File?         sequencing_reports           = sequencing_report.all_zip
        
        Array[String] data_tables_out              = select_first([data_tables.tables, []])
    }
}
