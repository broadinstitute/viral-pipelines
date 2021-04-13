version 1.0

import "sarscov2_illumina_full.wdl" as full_viral
import "../tasks/tasks_utils.wdl" as utils
# import "../tasks/tasks_terra.wdl" as terra

workflow workflow_launcher_sarscov2_illumina_full {
    meta {
        description: "Run full SARS-CoV-2 analysis workflow, capture outputs, write to new Terra data table."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        File          flowcell_tgz
        File          reference_fasta
        String        amplicon_bed_prefix

        Array[File]   biosample_attributes
        String        instrument_model
        String        sra_title

        # Int           min_genome_bases
        # Int           max_vadr_alerts

        String        workspace_name
        String        terra_project

        Array[File]+  samplesheets
        File          spikein_db

        String        account_name
        File          author_template_sbt
        String        spuid_namespace
    }
    String  flowcell_id = basename(basename(basename(basename(flowcell_tgz, ".gz"), ".zst"), ".tar"), ".tgz")

    call full_viral.sarscov2_illumina_full {
        input:
            flowcell_tgz = flowcell_tgz,
            reference_fasta = reference_fasta,
            amplicon_bed_prefix = amplicon_bed_prefix,
            biosample_attributes = biosample_attributes,
            instrument_model = instrument_model,
            sra_title = sra_title,
            workspace_name = workspace_name,
            terra_project = terra_project,
            samplesheets = samplesheets,
            spikein_db = spikein_db,
            account_name = account_name,
            author_template_sbt = author_template_sbt,
            spuid_namespace = spuid_namespace
    }

    # concatenate header and row tsv files
    # call utils.concatenate as flowcell_meta_tsv {
    #     input:
    #         infiles = [write_tsv([flowcell_tsv_header]), write_tsv([flowcell_tsv_row])],
    #         output_name = "flowcell_metadata-~{flowcell_id}.tsv"
    # }

     # upload concatenated .tsv file to flowcell table
    # call terra.upload_data_table_tsv as upload_flowcell_table {
    #     input:
    #         input_tsv = flowcell_meta_tsv.combined,
    #         terra_project = terra_project,
    #         workspace_name = workspace_name
    # }

    call gather_sarscov2_outputs {
        input:
            raw_reads_unaligned_bams     = sarscov2_illumina_full.raw_reads_unaligned_bams,
            cleaned_reads_unaligned_bams = sarscov2_illumina_full.cleaned_reads_unaligned_bams,
            cleaned_bams_tiny            = sarscov2_illumina_full.cleaned_bams_tiny,
            meta_by_filename_json        = sarscov2_illumina_full.meta_by_filename_json,
            read_counts_raw              = sarscov2_illumina_full.read_counts_raw,
            read_counts_depleted         = sarscov2_illumina_full.read_counts_depleted,
            sra_metadata                 = sarscov2_illumina_full.sra_metadata,
            cleaned_bam_uris             = sarscov2_illumina_full.cleaned_bam_uris,
            assemblies_fasta             = sarscov2_illumina_full.assemblies_fasta,
            passing_assemblies_fasta     = sarscov2_illumina_full.passing_assemblies_fasta,
            submittable_assemblies_fasta = sarscov2_illumina_full.submittable_assemblies_fasta,
            max_ntc_bases                = sarscov2_illumina_full.max_ntc_bases,
            demux_metrics                = sarscov2_illumina_full.demux_metrics,
            demux_commonBarcodes         = sarscov2_illumina_full.demux_commonBarcodes,
            demux_outlierBarcodes        = sarscov2_illumina_full.demux_outlierBarcodes,
            primer_trimmed_read_count    = sarscov2_illumina_full.primer_trimmed_read_count,
            primer_trimmed_read_percent  = sarscov2_illumina_full.primer_trimmed_read_percent,
            ivar_trim_stats_html         = sarscov2_illumina_full.ivar_trim_stats_html,
            ivar_trim_stats_png          = sarscov2_illumina_full.ivar_trim_stats_png,
            ivar_trim_stats_tsv          = sarscov2_illumina_full.ivar_trim_stats_tsv,
            multiqc_report_raw           = sarscov2_illumina_full.multiqc_report_raw,
            multiqc_report_cleaned       = sarscov2_illumina_full.multiqc_report_cleaned,
            spikein_counts               = sarscov2_illumina_full.spikein_counts,
            picard_metrics_wgs           = sarscov2_illumina_full.picard_metrics_wgs,
            assembly_stats_tsv           = sarscov2_illumina_full.assembly_stats_tsv,
            submission_zip               = sarscov2_illumina_full.submission_zip,
            submission_xml               = sarscov2_illumina_full.submission_xml,
            submit_ready                 = sarscov2_illumina_full.submit_ready,
            vadr_outputs                 = sarscov2_illumina_full.vadr_outputs,
            genbank_source_table         = sarscov2_illumina_full.genbank_source_table,
            gisaid_fasta                 = sarscov2_illumina_full.gisaid_fasta,
            gisaid_meta_tsv              = sarscov2_illumina_full.gisaid_meta_tsv,
            genbank_fasta                = sarscov2_illumina_full.genbank_fasta,
            nextmeta_tsv                 = sarscov2_illumina_full.nextmeta_tsv,
            nextclade_all_json           = sarscov2_illumina_full.nextclade_all_json,
            nextclade_auspice_json       = sarscov2_illumina_full.nextclade_auspice_json,
            assembled_ids                = sarscov2_illumina_full.assembled_ids,
            submittable_ids              = sarscov2_illumina_full.submittable_ids,
            failed_assembly_ids          = sarscov2_illumina_full.failed_assembly_ids,
            failed_annotation_ids        = sarscov2_illumina_full.failed_annotation_ids,
            num_read_files               = sarscov2_illumina_full.num_read_files,
            num_assembled                = sarscov2_illumina_full.num_assembled,
            num_failed_assembly          = sarscov2_illumina_full.num_failed_assembly,
            num_submittable              = sarscov2_illumina_full.num_submittable,
            num_failed_annotation        = sarscov2_illumina_full.num_failed_annotation,
            num_samples                  = sarscov2_illumina_full.num_samples,
            run_date                     = sarscov2_illumina_full.run_date,
            sequencing_reports           = sarscov2_illumina_full.sequencing_reports,
            data_tables_out              = sarscov2_illumina_full.data_tables_out
    }

    output {
        File flowcell_load_tsv = gather_sarscov2_outputs.flowcell_load_tsv
    }
}

task gather_sarscov2_outputs {
    input {
        Array[String]           raw_reads_unaligned_bams
        Array[String]           cleaned_reads_unaligned_bams
        Array[String]           cleaned_bams_tiny

        String                  meta_by_filename_json

        Array[Int]              read_counts_raw
        Array[Int]              read_counts_depleted

        String                  sra_metadata
        String                  cleaned_bam_uris

        Array[String]           assemblies_fasta
        Array[String]           passing_assemblies_fasta
        Array[String]           submittable_assemblies_fasta

        Int                     max_ntc_bases

        Array[String]           demux_metrics
        Array[String]           demux_commonBarcodes
        Array[String]           demux_outlierBarcodes

        Array[Int]              primer_trimmed_read_count
        Array[Float]            primer_trimmed_read_percent
        String                  ivar_trim_stats_html
        String                  ivar_trim_stats_png
        String                  ivar_trim_stats_tsv
        String                  multiqc_report_raw
        String                  multiqc_report_cleaned
        String                  spikein_counts
        String                  picard_metrics_wgs

        String                  assembly_stats_tsv

        String                  submission_zip
        String                  submission_xml
        String                  submit_ready
        Array[String]           vadr_outputs
        String                  genbank_source_table

        String                  gisaid_fasta
        String                  gisaid_meta_tsv

        String                  genbank_fasta
        String                  nextmeta_tsv

        String                  nextclade_all_json
        String                  nextclade_auspice_json

        Array[String]           assembled_ids
        Array[String]           submittable_ids
        Array[String]           failed_assembly_ids
        Array[String]           failed_annotation_ids
        Int                     num_read_files
        Int                     num_assembled
        Int                     num_failed_assembly
        Int                     num_submittable
        Int                     num_failed_annotation
        Int                     num_samples

        String                  run_date

        String?                 sequencing_reports

        Array[String]           data_tables_out
    }

    command <<<
        # create header line in final output load file
        echo -e "entity:flowcell_id\tassembled_ids\tassemblies_fasta\tassembly_stats_tsv\tauthors_sbt\t \
                cleaned_bam_uris\tcleaned_reads_unaligned_bams\tdata_table_status\t \
                demux_commonBarcodes\tdemux_metrics\tdemux_outlierBarcodes\t \
                failed_annotation_ids\tfailed_assembly_ids\tflowcell_tgz\tgenbank_fasta\t \
                genbank_source_table\tgisaid_fasta\tgisaid_meta_tsv\tinstrument_model\t \
                ivar_trim_stats_html\tivar_trim_stats_png\tivar_trim_stats_tsv\tmax_ntc_bases\t \
                meta_by_filename_json\tmultiqc_report_cleaned\tmultiqc_report_raw\t \
                nextclade_all_json\tnextclade_auspice_json\tnextmeta_tsv\tnum_assembled\t \
                num_failed_annotation\tnum_failed_assembly\tnum_read_files\tnum_samples\t \
                num_submittable\tpassing_assemblies_fasta\tpicard_metrics_wgs\t \
                primer_trimmed_read_count\tprimer_trimmed_read_percent\traw_reads_unaligned_bams\t \
                read_counts_depleted\tread_counts_raw\trun_date\tsamplesheets\tsequencing_reports\t \
                spikein_counts\tsra_metadata\tsubmission_xml\tsubmission_zip\tsubmit_ready\t \
                submittable_assemblies_fasta\tsubmittable_ids\ttitle\tvadr_outputs" > flowcell_load_table.tsv

        echo -e "my_test_unique_identifier_snapshot\t \
                ~{sep="," raw_reads_unaligned_bams}\t \
                ~{sep="," cleaned_reads_unaligned_bams}\t \
                ~{sep="," cleaned_bams_tiny}\t \
                ~{meta_by_filename_json}\t \
                ~{sep="," read_counts_raw}\t \
                ~{sep="," read_counts_depleted}\t \
                ~{sra_metadata}\t \
                ~{cleaned_bam_uris}\t \
                ~{sep="," assemblies_fasta}\t \
                ~{sep="," passing_assemblies_fasta}\t \
                ~{sep="," submittable_assemblies_fasta}\t \
                ~{max_ntc_bases}\t \
                ~{sep="," demux_metrics}\t \
                ~{sep="," demux_commonBarcodes}\t \
                ~{sep="," demux_outlierBarcodes}\t \
                ~{sep="," primer_trimmed_read_count}\t \
                ~{sep="," primer_trimmed_read_percent}\t \
                ~{ivar_trim_stats_html}\t \
                ~{ivar_trim_stats_png}\t \
                ~{ivar_trim_stats_tsv}\t \
                ~{multiqc_report_raw}\t \
                ~{multiqc_report_cleaned}\t \
                ~{spikein_counts}\t \
                ~{picard_metrics_wgs}\t \
                ~{assembly_stats_tsv}\t \
                ~{submission_zip}\t \
                ~{submission_xml}\t \
                ~{submit_ready}\t \
                ~{sep="," vadr_outputs}\t \
                ~{genbank_source_table}\t \
                ~{gisaid_fasta}\t \
                ~{gisaid_meta_tsv}\t \
                ~{genbank_fasta}\t \
                ~{nextmeta_tsv}\t \
                ~{nextclade_all_json}\t \
                ~{nextclade_auspice_json}\t \
                ~{sep="," assembled_ids}\t \
                ~{sep="," submittable_ids}\t \
                ~{sep="," failed_assembly_ids}\t \
                ~{sep="," failed_annotation_ids}\t \
                ~{num_read_files}\t \
                ~{num_assembled}\t \
                ~{num_failed_assembly}\t \
                ~{num_submittable}\t \
                ~{num_failed_annotation}\t \
                ~{num_samples}\t \
                ~{run_date}\t \
                ~{sequencing_reports}\t \
                ~{sep=","data_tables_out}" >> flowcell_load_table.tsv
    >>>

    runtime {
        docker: "quay.io/broadinstitute/viral-core:2.1.19"

    }

    output {
        File flowcell_load_tsv = "flowcell_load_table.tsv"
    }
}
