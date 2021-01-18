version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_demux.wdl" as demux

import "assemble_refbased.wdl"
import "sarscov2_lineages.wdl"

workflow sarscov2_illumina_full {
    meta {
        description: "Full SARS-CoV-2 analysis workflow starting from raw Illumina flowcell (tar.gz) and metadata and performing assembly, spike-in analysis, qc, lineage assignment, and packaging for data release."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
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
        File          flowcell_tgz
        Array[File]+  samplesheets  ## must be in lane order!

        File          reference_fasta
        String        amplicon_bed_prefix

        File          biosample_attributes
        File?         sample_rename_map

        Int           taxid = 2697049
        Int           min_genome_bases = 20000
        Int           min_reads_per_bam = 100
        String        gisaid_prefix = 'hCoV-19/'

        File          spikein_db
        Array[File]?  bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]?  blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]?  bwaDbs
    }

    #### demux each lane (rename samples if requested)
    scatter(lane_sheet in zip(range(length(samplesheets)), samplesheets)) {
        call demux.samplesheet_rename_ids {
            input:
                old_sheet = lane_sheet.right,
                rename_map = sample_rename_map
        }
        call demux.illumina_demux {
            input:
                flowcell_tgz = flowcell_tgz,
                lane = lane_sheet.left + 1,
                samplesheet = samplesheet_rename_ids.new_sheet
        }
        call demux.map_map_setdefault as meta_default_sample {
            input:
                map_map_json = illumina_demux.meta_by_sample_json,
                sub_keys = ["amplicon_set", "control"]
        }
        call demux.map_map_setdefault as meta_default_filename {
            input:
                map_map_json = illumina_demux.meta_by_filename_json,
                sub_keys = ["spike_in"]
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

        # TO DO: flag all libraries where highest spike-in is not what was expected in extended samplesheet
    }

    #### SRA submission prep
    call ncbi.sra_meta_prep {
        input:
            cleaned_bam_filepaths = select_all(cleaned_bam_passing),
            biosample_map = biosample_attributes,
            library_metadata = samplesheet_rename_ids.new_sheet,
            out_name = "sra_metadata-~{basename(flowcell_tgz, '.tar.gz')}.tsv",
            platform = "ILLUMINA"
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

    ### gather data by biosample
    call read_utils.group_bams_by_sample {
        input:
            bam_filepaths = select_all(cleaned_bam_passing)
    }

    ### assembly and analyses per biosample
    scatter(name_reads in zip(group_bams_by_sample.sample_names, group_bams_by_sample.grouped_bam_filepaths)) {
        Boolean ampseq = (meta_sample.merged[name_reads.left]["amplicon_set"] != "")
        String orig_name = meta_sample.merged[name_reads.left]["sample_original"]

        # assemble genome
        if (ampseq) {
            String trim_coords_bed = amplicon_bed_prefix + meta_sample.merged[name_reads.left]["amplicon_set"] + ".bed"
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
        if (meta_sample.merged[name_reads.left]["control"] == 'NTC') {
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
            Array[String] assembly_cmt = [orig_name, "Broad viral-ngs v. " + illumina_demux.viralngs_version[0], assemble_refbased.assembly_mean_coverage]

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

        Map[String,String?] assembly_stats = {
            'sample_orig': orig_name,
            'sample': name_reads.left,
            'amplicon_set': meta_sample.merged[name_reads.left]["amplicon_set"],
            'assembly_mean_coverage': assemble_refbased.assembly_mean_coverage,
            'nextclade_clade':   sarscov2_lineages.nextclade_clade,
            'nextclade_aa_subs': sarscov2_lineages.nextclade_aa_subs,
            'nextclade_aa_dels': sarscov2_lineages.nextclade_aa_dels,
            'pango_lineage':     sarscov2_lineages.pango_lineage
        }
        Map[String,File?] assembly_files = {
            'assembly_fasta':           assemble_refbased.assembly_fasta,
            'coverage_plot':            assemble_refbased.align_to_ref_merged_coverage_plot,
            'aligned_bam':              assemble_refbased.align_to_ref_merged_aligned_trimmed_only_bam,
            'replicate_discordant_vcf': assemble_refbased.replicate_discordant_vcf,
            'nextclade_tsv': sarscov2_lineages.nextclade_tsv,
            'pangolin_csv':  sarscov2_lineages.pangolin_csv,
            'vadr_tgz': vadr.outputs_tgz
        }
        Map[String,Int?] assembly_metrics = {
            'assembly_length_unambiguous': assemble_refbased.assembly_length_unambiguous,
            'dist_to_ref_snps':            assemble_refbased.dist_to_ref_snps,
            'dist_to_ref_indels':          assemble_refbased.dist_to_ref_indels,
            'replicate_concordant_sites':  assemble_refbased.replicate_concordant_sites,
            'replicate_discordant_snps':   assemble_refbased.replicate_discordant_snps,
            'replicate_discordant_indels': assemble_refbased.replicate_discordant_indels,
            'num_read_groups':             assemble_refbased.num_read_groups,
            'num_libraries':               assemble_refbased.num_libraries,
            'vadr_num_alerts': vadr.num_alerts
        }

    }

    # TO DO: filter out genomes from submission that are less than ntc_bases.max
    call read_utils.max as ntc {
      input:
        list = select_all(ntc_bases)
    }

    ### prep genbank submission
    call nextstrain.concatenate as submit_genomes {
      input:
        infiles = select_all(submittable_genomes),
        output_name = "assemblies.fasta"
    }
    call ncbi.biosample_to_genbank {
      input:
        biosample_attributes = biosample_attributes,
        num_segments = 1,
        taxid = taxid,
        filter_to_ids = write_lines(select_all(submittable_id))
    }
    call ncbi.structured_comments {
      input:
        assembly_stats_tsv = write_tsv(flatten([[['SeqID','Assembly Method','Coverage']],select_all(assembly_cmt)])),
        filter_to_ids = write_lines(select_all(submittable_id))
    }
    call ncbi.package_genbank_ftp_submission {
      input:
        sequences_fasta = submit_genomes.combined,
        source_modifier_table = biosample_to_genbank.genbank_source_modifier_table,
        structured_comment_table = structured_comments.structured_comment_table
    }

    ### prep gisaid submission
    call ncbi.prefix_fasta_header as prefix_gisaid {
      input:
        genome_fasta = submit_genomes.combined,
        prefix = gisaid_prefix
    }
    call ncbi.gisaid_meta_prep {
      input:
        source_modifier_table = biosample_to_genbank.genbank_source_modifier_table,
        structured_comments = structured_comments.structured_comment_table,
        out_name = "gisaid_meta.tsv"
    }

    output {
        Array[File] raw_reads_unaligned_bams     = flatten(illumina_demux.raw_reads_unaligned_bams)
        Array[File] cleaned_reads_unaligned_bams = select_all(cleaned_bam_passing)
        Array[File] cleaned_bams_tiny = select_all(empty_bam)

        File meta_by_filename = meta_filename.merged_json

        Array[Int]  read_counts_raw = deplete.depletion_read_count_pre
        Array[Int]  read_counts_depleted = deplete.depletion_read_count_post

        File        sra_metadata          = sra_meta_prep.sra_metadata

        Array[File] assemblies_fasta = assemble_refbased.assembly_fasta
        Array[File] passing_assemblies_fasta = select_all(passing_assemblies)
        Array[File] submittable_assemblies_fasta = select_all(submittable_genomes)

        Int         max_ntc_bases = ntc.max

        Array[File] demux_metrics            = illumina_demux.metrics
        Array[File] demux_commonBarcodes     = illumina_demux.commonBarcodes
        Array[File] demux_outlierBarcodes    = illumina_demux.outlierBarcodes

        File        multiqc_report_raw     = multiqc_raw.multiqc_report
        File        multiqc_report_cleaned = multiqc_cleaned.multiqc_report
        File        spikein_counts         = spike_summary.count_summary

        Array[Map[String,String?]] per_assembly_stats = assembly_stats
        Array[Map[String,File?]]   per_assembly_files = assembly_files
        Array[Map[String,Int?]]    per_assembly_metrics = assembly_metrics

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
        Int           num_read_files = length(select_all(cleaned_bam_passing))
        Int           num_assembled = length(select_all(passing_assemblies))
        Int           num_failed_assembly = length(select_all(failed_assembly_id))
        Int           num_submittable = length(select_all(submittable_id))
        Int           num_failed_annotation = length(select_all(failed_annotation_id))
        Int           num_samples = length(group_bams_by_sample.sample_names)
    }
}
