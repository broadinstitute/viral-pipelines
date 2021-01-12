version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_sarscov2.wdl" as sarscov2

import "assemble_refbased.wdl"
import "sarscov2_lineages.wdl"
import "sarscov2_genbank.wdl"
import "../tasks/tasks_demux.wdl" as demux

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

        authors_sbt: {
          description: "A genbank submission template file (SBT) with the author list, created at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/",
          patterns: ["*.sbt"]
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

        File          authors_sbt
        File          biosample_attributes
        File?         fasta_rename_map

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
        call assemble_refbased.assemble_refbased {
            input:
                reads_unmapped_bams = [name_reads.right],
                reference_fasta = reference_fasta,
                sample_name = name_reads.left
                # lookup skip_mark_dupes and trim_coords_bed from metadata
        }

        if ( assemble_refbased.assembly_length_unambiguous >= min_genome_bases) {
            File passing_assemblies = assemble_refbased.assembly_fasta
            call sarscov2_lineages.sarscov2_lineages {
                input:
                    genome_fasta = assemble_refbased.assembly_fasta
            }
        }
    }

    ### prep genbank submission
    call sarscov2_genbank.sarscov2_genbank {
        input:
            assemblies_fasta = assemble_refbased.assembly_fasta,
            taxid = taxid,
            min_genome_bases = min_genome_bases
    }

    output {
        Array[File] raw_reads_unaligned_bams     = illumina_demux.raw_reads_unaligned_bams
        Array[File] cleaned_reads_unaligned_bams = deplete.cleaned_bam

        Array[Int]  read_counts_raw = deplete.depletion_read_count_pre
        Array[Int]  read_counts_depleted = deplete.depletion_read_count_post

        File        sra_metadata          = sra_meta_prep.sra_metadata

        Array[File] assemblies_fasta = assemble_refbased.assembly_fasta

        File        demux_metrics            = illumina_demux.metrics
        File        demux_commonBarcodes     = illumina_demux.commonBarcodes
        File        demux_outlierBarcodes    = illumina_demux.outlierBarcodes

        File        multiqc_report_raw     = multiqc_raw.multiqc_report
        File        multiqc_report_cleaned = multiqc_cleaned.multiqc_report
        File        spikein_counts         = spike_summary.count_summary
    }
}
