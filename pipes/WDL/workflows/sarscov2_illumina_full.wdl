version 1.0

import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_reports.wdl" as reports

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
        samplesheets: {
            description: "Custom formatted 'extended' format tsv samplesheets that will override any SampleSheet.csv in the illumina BCL directory. Must supply one file per lane of the flowcell, and must provide them in lane order. Required tsv column headings are: sample, library_id_per_sample, barcode_1, barcode_2 (if paired reads, omit if single-end), library_strategy, library_source, library_selection, design_description. 'sample' must correspond to a biological sample. 'sample' x 'library_id_per_sample' must be unique within a samplesheet and correspond to independent libraries from the same original sample. barcode_1 and barcode_2 must correspond to the actual index sequence. Remaining columns must follow strict ontology: see 3rd tab of https://www.ncbi.nlm.nih.gov/core/assets/sra/files/SRA_metadata_acc_example.xlsx for controlled vocabulary and term definitions.",
            patterns: ["*.txt", "*.tsv"]
        }

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

    }

    input {
        File          reference_fasta
        String        amplicon_bed_prefix

        Array[File]   biosample_attributes
        String        instrument_model
        String        sra_title

        Int           min_genome_bases = 15000
    }
    Int     taxid = 2697049
    String  gisaid_prefix = 'hCoV-19/'

    # merge biosample attributes tables
    call reports.tsv_join as biosample_merge {
        input:
            input_tsvs = biosample_attributes,
            id_col = 'accession',
            out_basename = "biosample_attributes-merged"
    }

    ### demux, deplete, SRA submission prep, fastqc/multiqc
    call demux_deplete.demux_deplete {
        input:
            biosample_map = biosample_merge.out_tsv,
            instrument_model = instrument_model,
            sra_title = sra_title
    }

    ### gather data by biosample
    call read_utils.group_bams_by_sample {
        input:
            bam_filepaths = demux_deplete.cleaned_reads_unaligned_bams
    }

    ### assembly and analyses per biosample
    scatter(name_reads in zip(group_bams_by_sample.sample_names, group_bams_by_sample.grouped_bam_filepaths)) {
        Boolean ampseq = (demux_deplete.meta_by_sample[name_reads.left]["amplicon_set"] != "")
        String orig_name = demux_deplete.meta_by_sample[name_reads.left]["sample_original"]

        # assemble genome
        if (ampseq) {
            String trim_coords_bed = amplicon_bed_prefix + demux_deplete.meta_by_sample[name_reads.left]["amplicon_set"] + ".bed"
        }
        call assemble_refbased.assemble_refbased {
            input:
                reads_unmapped_bams = name_reads.right,
                reference_fasta = reference_fasta,
                sample_name = name_reads.left,
                aligner = "minimap2",
                skip_mark_dupes = ampseq,
                trim_coords_bed = trim_coords_bed,
                min_coverage = if ampseq then 20 else 3
        }

        # log controls
        if (demux_deplete.meta_by_sample[name_reads.left]["control"] == 'NTC') {
            Int ntc_bases = assemble_refbased.assembly_length_unambiguous
        }

        # for genomes that somewhat assemble
        if (assemble_refbased.assembly_length_unambiguous >= min_genome_bases) {
            call ncbi.rename_fasta_header {
              input:
                genome_fasta = assemble_refbased.assembly_fasta,
                new_name = orig_name
            }

            File passing_assemblies = rename_fasta_header.renamed_fasta
            String passing_assembly_ids = orig_name
            Array[String] assembly_cmt = [orig_name, "Broad viral-ngs v. " + demux_deplete.demux_viral_core_version, assemble_refbased.assembly_mean_coverage, instrument_model]

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
            if (vadr.num_alerts==0) {
              File submittable_genomes = passing_assemblies
              String submittable_id = orig_name
            }
            if (vadr.num_alerts>0) {
              String failed_annotation_id = orig_name
            }
        }
        if (assemble_refbased.assembly_length_unambiguous < min_genome_bases) {
            String failed_assembly_id = orig_name
        }

        Array[String] assembly_tsv_row = [
            orig_name,
            name_reads.left,
            demux_deplete.meta_by_sample[name_reads.left]["amplicon_set"],
            assemble_refbased.assembly_mean_coverage,
            assemble_refbased.assembly_length_unambiguous,
            select_first([sarscov2_lineages.nextclade_clade, ""]),
            select_first([sarscov2_lineages.nextclade_aa_subs, ""]),
            select_first([sarscov2_lineages.nextclade_aa_dels, ""]),
            select_first([sarscov2_lineages.pango_lineage, ""]),
            assemble_refbased.dist_to_ref_snps,
            assemble_refbased.dist_to_ref_indels,
            select_first([vadr.num_alerts, ""]),
            assemble_refbased.assembly_fasta,
            assemble_refbased.align_to_ref_merged_coverage_plot,
            assemble_refbased.align_to_ref_merged_aligned_trimmed_only_bam,
            assemble_refbased.replicate_discordant_vcf,
            select_first([sarscov2_lineages.nextclade_tsv, ""]),
            select_first([sarscov2_lineages.pangolin_csv, ""]),
            select_first([vadr.outputs_tgz, ""]),
            assemble_refbased.replicate_concordant_sites,
            assemble_refbased.replicate_discordant_snps,
            assemble_refbased.replicate_discordant_indels,
            assemble_refbased.num_read_groups,
            assemble_refbased.num_libraries,
        ]
    }
    Array[String] assembly_tsv_header = [
        'sample', 'sample_sanitized', 'amplicon_set', 'assembly_mean_coverage', 'assembly_length_unambiguous',
        'nextclade_clade', 'nextclade_aa_subs', 'nextclade_aa_dels', 'pango_lineage',
        'dist_to_ref_snps', 'dist_to_ref_indels', 'vadr_num_alerts',
        'assembly_fasta', 'coverage_plot', 'aligned_bam', 'replicate_discordant_vcf',
        'nextclade_tsv', 'pangolin_csv', 'vadr_tgz',
        'replicate_concordant_sites', 'replicate_discordant_snps', 'replicate_discordant_indels', 'num_read_groups', 'num_libraries',
        ]

    call nextstrain.concatenate as assembly_meta_tsv {
      input:
        infiles = [write_tsv([assembly_tsv_header]), write_tsv(assembly_tsv_row)],
        output_name = "assembly_metadata.tsv"
    }


    # TO DO: filter out genomes from submission that are less than ntc_max.out
    call read_utils.max as ntc_max {
      input:
        list = select_all(ntc_bases)
    }

    ### prep genbank submission
    call ncbi.biosample_to_genbank {
      input:
        biosample_attributes = biosample_merge.out_tsv,
        num_segments = 1,
        taxid = taxid,
        filter_to_ids = write_lines(select_all(submittable_id))
    }
    call ncbi.structured_comments {
      input:
        assembly_stats_tsv = write_tsv(flatten([[['SeqID','Assembly Method','Coverage','Sequencing Technology']],select_all(assembly_cmt)])),
        filter_to_ids = biosample_to_genbank.sample_ids
    }
    call nextstrain.concatenate as passing_genomes {
      input:
        infiles = select_all(submittable_genomes),
        output_name = "assemblies.fasta"
    }
    call nextstrain.filter_sequences_to_list as submit_genomes {
      input:
        sequences = passing_genomes.combined,
        keep_list = [biosample_to_genbank.sample_ids]
    }
    call ncbi.package_genbank_ftp_submission {
      input:
        sequences_fasta = submit_genomes.filtered_fasta,
        source_modifier_table = biosample_to_genbank.genbank_source_modifier_table,
        structured_comment_table = structured_comments.structured_comment_table
    }

    ### prep gisaid submission
    call ncbi.prefix_fasta_header as prefix_gisaid {
      input:
        genome_fasta = submit_genomes.filtered_fasta,
        prefix = gisaid_prefix
    }
    call ncbi.gisaid_meta_prep {
      input:
        source_modifier_table = biosample_to_genbank.genbank_source_modifier_table,
        structured_comments = structured_comments.structured_comment_table,
        out_name = "gisaid_meta.tsv"
    }

    output {
        Array[File] raw_reads_unaligned_bams     = demux_deplete.raw_reads_unaligned_bams
        Array[File] cleaned_reads_unaligned_bams = demux_deplete.cleaned_reads_unaligned_bams
        Array[File] cleaned_bams_tiny            = demux_deplete.cleaned_bams_tiny

        File meta_by_filename_json = demux_deplete.meta_by_filename_json

        Array[Int]  read_counts_raw       = demux_deplete.read_counts_raw
        Array[Int]  read_counts_depleted  = demux_deplete.read_counts_depleted

        File        sra_metadata          = select_first([demux_deplete.sra_metadata])

        Array[File] assemblies_fasta = assemble_refbased.assembly_fasta
        Array[File] passing_assemblies_fasta = select_all(passing_assemblies)
        Array[File] submittable_assemblies_fasta = select_all(submittable_genomes)

        Int         max_ntc_bases = ntc_max.out

        Array[File] demux_metrics            = demux_deplete.demux_metrics
        Array[File] demux_commonBarcodes     = demux_deplete.demux_commonBarcodes
        Array[File] demux_outlierBarcodes    = demux_deplete.demux_outlierBarcodes

        File        multiqc_report_raw     = demux_deplete.multiqc_report_raw
        File        multiqc_report_cleaned = demux_deplete.multiqc_report_cleaned
        File        spikein_counts         = demux_deplete.spikein_counts

        File assembly_stats_tsv = assembly_meta_tsv.combined

        File submission_zip = package_genbank_ftp_submission.submission_zip
        File submission_xml = package_genbank_ftp_submission.submission_xml
        File submit_ready   = package_genbank_ftp_submission.submit_ready
        Array[File] vadr_outputs = select_all(vadr.outputs_tgz)
        File genbank_source_table = biosample_to_genbank.genbank_source_modifier_table

        File gisaid_fasta = prefix_gisaid.renamed_fasta
        File gisaid_meta_tsv = gisaid_meta_prep.meta_tsv

        Array[String] assembled_ids = select_all(passing_assembly_ids)
        Array[String] submittable_ids = select_all(submittable_id)
        Array[String] failed_assembly_ids = select_all(failed_assembly_id)
        Array[String] failed_annotation_ids = select_all(failed_annotation_id)
        Int           num_read_files = length(demux_deplete.cleaned_reads_unaligned_bams)
        Int           num_assembled = length(select_all(passing_assemblies))
        Int           num_failed_assembly = length(select_all(failed_assembly_id))
        Int           num_submittable = length(select_all(submittable_id))
        Int           num_failed_annotation = length(select_all(failed_annotation_id))
        Int           num_samples = length(group_bams_by_sample.sample_names)
    }
}
