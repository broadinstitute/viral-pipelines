version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_sarscov2.wdl" as sarscov2
import "../tasks/tasks_demux.wdl" as demux

import "assemble_refbased.wdl"

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
        ampseq_trim_coords_bed: {
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

        # TO DO: flag all libraries where highest spike-in is not what was expected in extended samplesheet
    }

    #### SRA submission prep
    call ncbi.sra_meta_prep {
        input:
            cleaned_bam_filepaths = deplete.cleaned_bam,
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
            bam_filepaths = deplete.cleaned_bam
    }
    call read_utils.get_sample_meta {
        input:
            sanitized_sample_names = group_bams_by_sample.sample_names,
            samplesheets_extended = samplesheet_rename_ids.new_sheet
    }

    ### assembly and analyses per biosample
    scatter(name_reads in zip(group_bams_by_sample.sample_names, group_bams_by_sample.grouped_bam_filepaths)) {
        # assemble genome
        if (get_sample_meta.amplicon_set[name_reads.left] != "") {
            String trim_coords_bed = amplicon_bed_prefix + get_sample_meta.amplicon_set[name_reads.left] + ".bed"
        }
        call assemble_refbased.assemble_refbased {
            input:
                reads_unmapped_bams = name_reads.right,
                reference_fasta = reference_fasta,
                sample_name = name_reads.left,
                aligner = "minimap2",
                skip_mark_dupes = (get_sample_meta.amplicon_set[name_reads.left] != ""),
                trim_coords_bed = trim_coords_bed,
                min_coverage = if (get_sample_meta.amplicon_set[name_reads.left] != "") then 20 else 3
        }

        # log controls
        if (get_sample_meta.control[name_reads.left] == 'NTC') {
            Int ntc_bases = assemble_refbased.assembly_length_unambiguous
        }

        # for genomes that somewhat assemble
        if (assemble_refbased.assembly_length_unambiguous >= min_genome_bases) {
            call ncbi.rename_fasta_header {
              input:
                genome_fasta = assemble_refbased.assembly_fasta,
                new_name = get_sample_meta.original_names[name_reads.left]
            }

            File passing_assemblies = rename_fasta_header.renamed_fasta
            String passing_assembly_ids = get_sample_meta.original_names[name_reads.left]
            Array[String] assembly_meta = [get_sample_meta.original_names[name_reads.left], assemble_refbased.assembly_mean_coverage]

            # lineage assignment
            call sarscov2.nextclade_one_sample {
                input:
                    genome_fasta = passing_assemblies
            }
            call sarscov2.pangolin_one_sample {
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
              String submittable_id = name_reads.left
            }
            if (vadr.num_alerts>0) {
              String failed_annotation_id = name_reads.left
            }
        }
        if (assemble_refbased.assembly_length_unambiguous < min_genome_bases) {
            String failed_assembly_id = name_reads.left
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
        assembly_stats_tsv = write_tsv(flatten([[['SeqID','Coverage']],select_all(assembly_meta)])),
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

    output {
        Array[File] raw_reads_unaligned_bams     = flatten(illumina_demux.raw_reads_unaligned_bams)
        Array[File] cleaned_reads_unaligned_bams = deplete.cleaned_bam

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

        # TO DO: bundle outputs into structs or some meaningful thing
        #String nextclade_clade = nextclade_one_sample.nextclade_clade
        #File   nextclade_tsv   = nextclade_one_sample.nextclade_tsv
        #String nextclade_aa_subs = nextclade_one_sample.aa_subs_csv
        #String nextclade_aa_dels = nextclade_one_sample.aa_dels_csv
        #String pangolin_clade  = pangolin_one_sample.pangolin_clade
        #File   pangolin_csv    = pangolin_one_sample.pangolin_csv

        File submission_zip = package_genbank_ftp_submission.submission_zip
        File submission_xml = package_genbank_ftp_submission.submission_xml
        File submit_ready   = package_genbank_ftp_submission.submit_ready
        Array[File] vadr_outputs = select_all(vadr.outputs_tgz)
        File genbank_source_table = biosample_to_genbank.genbank_source_modifier_table

        File gisaid_fasta = prefix_gisaid.renamed_fasta

        Array[String] assembled_ids = select_all(passing_assembly_ids)
        Array[String] submittable_ids = select_all(submittable_id)
        Array[String] failed_assembly_ids = select_all(failed_assembly_id)
        Array[String] failed_annotation_ids = select_all(failed_annotation_id)
        Int           num_assembled = length(select_all(passing_assemblies))
        Int           num_failed_assembly = length(select_all(failed_assembly_id))
        Int           num_submittable = length(select_all(submittable_id))
        Int           num_failed_annotation = length(select_all(failed_annotation_id))
        Int           num_samples = length(group_bams_by_sample.sample_names)
    }
}
