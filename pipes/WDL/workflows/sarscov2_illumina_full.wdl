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
        reference_fasta: {
            description: "Reference genome to align reads to.",
            patterns: ["*.fasta"]
        }
        ampseq_trim_coords_bed: {
            description: "amplicon primers to trim in reference coordinate space (0-based BED format)",
            patterns: ["*.bed"]
        }

        biosample_attributes: {
          description: "A post-submission attributes file from NCBI BioSample, which is available at https://submit.ncbi.nlm.nih.gov/subs/ and clicking on 'Download attributes file with BioSample accessions'.",
          patterns: ["*.txt", "*.tsv"]
        }

    }

    input {
        File         flowcell_tgz
        Array[File]+ samplesheets  ## must be in lane order!

        File          reference_fasta
        File          ampseq_trim_coords_bed

        File          biosample_attributes
        File?         rename_map

        Int           taxid = 2697049
        Int           min_genome_bases = 20000

        File          spikein_db
        File          trim_clip_db
        Array[File]?  bmtaggerDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]?  blastDbs  # .tar.gz, .tgz, .tar.bz2, .tar.lz4, .fasta, or .fasta.gz
        Array[File]?  bwaDbs
    }

    #### demux each lane
    scatter(lane_sheet in zip(range(length(samplesheets)), samplesheets)) {
        call demux.illumina_demux as illumina_demux {
            input:
                flowcell_tgz = flowcell_tgz,
                lane = lane_sheet.left + 1,
                samplesheet = lane_sheet.right
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
    }

    #### SRA submission prep
    call ncbi.sra_meta_prep {
        input:
            cleaned_bam_filepaths = deplete.cleaned_bam,
            out_name = "sra_metadata-~{basename(flowcell_tgz, '.tar.gz')}.tsv"
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

    ### assembly and analyses per biosample
    call read_utils.group_bams_by_sample {
        input:
            bam_filepaths = deplete.cleaned_bam
    }

    scatter(name_reads in zip(group_bams_by_sample.sample_names, group_bams_by_sample.grouped_bam_filepaths)) {
        # assemble genome
        call assemble_refbased.assemble_refbased {
            input:
                reads_unmapped_bams = name_reads.right,
                reference_fasta = reference_fasta,
                sample_name = name_reads.left
                # TO DO: lookup skip_mark_dupes and trim_coords_bed from metadata
        }

        # for genomes that somewhat assemble
        if (assemble_refbased.assembly_length_unambiguous >= min_genome_bases) {
            File passing_assemblies = assemble_refbased.assembly_fasta
            String passing_assembly_ids = name_reads.left
            Array[String] assembly_meta = [name_reads.left, assemble_refbased.assembly_mean_coverage]

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

    # TO DO: add some error checks / filtration if NTCs assemble above some length or above other genomes

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

    output {
        Array[File] raw_reads_unaligned_bams     = flatten(illumina_demux.raw_reads_unaligned_bams)
        Array[File] cleaned_reads_unaligned_bams = deplete.cleaned_bam

        Array[Int]  read_counts_raw = deplete.depletion_read_count_pre
        Array[Int]  read_counts_depleted = deplete.depletion_read_count_post

        File        sra_metadata          = sra_meta_prep.sra_metadata

        Array[File] assemblies_fasta = assemble_refbased.assembly_fasta
        Array[File] passing_assemblies_fasta = select_all(passing_assemblies)
        Array[File] submittable_assemblies_fasta = select_all(submittable_genomes)

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
