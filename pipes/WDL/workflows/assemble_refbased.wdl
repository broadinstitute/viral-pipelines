version 1.0

import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_read_utils.wdl" as read_utils

workflow assemble_refbased {

    meta {
        description: "Reference-based microbial consensus calling. Aligns NGS reads to a singular reference genome, calls a new consensus sequence, and emits: new assembly, reads aligned to provided reference, reads aligned to new assembly, various figures of merit, plots, and QC metrics. The user may provide unaligned reads spread across multiple input files and this workflow will parallelize alignment per input file before merging results prior to consensus calling."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    parameter_meta {
        sample_name: {
            description: "Base name of output files. The 'SM' field in BAM read group headers are also rewritten to this value. Avoid spaces and other filename-unfriendly characters.",
            category: "common"
        }
        reads_unmapped_bams: {
            description: "Unaligned reads in BAM format",
            patterns: ["*.bam"]
        }
        reference_fasta: {
            description: "Reference genome to align reads to.",
            patterns: ["*.fasta"]
        }
        aligner: {
            description: "Read aligner software to use. Options: novoalign, bwa, minimap2. Minimap2 can automatically handle Illumina, PacBio, or Oxford Nanopore reads as long as the 'PL' field in the BAM read group header is set properly (novoalign and bwa are Illumina-only)."
        }
        novocraft_license: {
            description: "The default Novoalign short read aligner is a commercially licensed software that is available in a much slower, single-threaded version for free. If you have a paid license file, provide it here to run in multi-threaded mode. If this is omitted, it will run in single-threaded mode.",
            patterns: ["*.lic"]
        }
        skip_mark_dupes: {
            description: "skip Picard MarkDuplicates step after alignment. This is recommended to be set to true for PCR amplicon based data. (Default: false)"
        }
        trim_coords_bed: {
            description: "optional primers to trim in reference coordinate space (0-based BED format)",
            patterns: ["*.bed"],
            category: "common"
        }


        assembly_fasta: { description: "The new assembly / consensus sequence for this sample" }
        align_to_ref_variants_vcf_gz: { description: "All variants in the input reads against the original reference genome. This VCF file is used to create the assembly_fasta" }
        assembly_length: { description: "The length of the sequence described in assembly_fasta, inclusive of any uncovered regions denoted by Ns" }
        assembly_length_unambiguous: { description: "The number of called consensus bases in assembly_fasta (excludes regions of the genome that lack read coverage)" }
    }

    input {
        String          sample_name = basename(reads_unmapped_bams[0], '.bam')
        Array[File]+    reads_unmapped_bams
        File            reference_fasta

        String          aligner="minimap2"
        File?           novocraft_license
        Boolean?        skip_mark_dupes=false
        File?           trim_coords_bed
    }

    Map[String,String] align_to_ref_options = {
                            "novoalign": "-r Random -l 40 -g 40 -x 20 -t 501 -k",
                            "bwa": "-k 12 -B 1",
                            "minimap2": ""
                            }
    Map[String,String] align_to_self_options = {
                            "novoalign": "-r Random -l 40 -g 40 -x 20 -t 100",
                            "bwa": "",
                            "minimap2": ""
                            }

    scatter(reads_unmapped_bam in reads_unmapped_bams) {
        call assembly.align_reads as align_to_ref {
            input:
                reference_fasta    = reference_fasta,
                reads_unmapped_bam = reads_unmapped_bam,
                novocraft_license  = novocraft_license,
                skip_mark_dupes    = skip_mark_dupes,
                aligner            = aligner,
                aligner_options    = align_to_ref_options[aligner]
        }
        call assembly.ivar_trim {
            input:
                aligned_bam = align_to_ref.aligned_only_reads_bam,
                trim_coords_bed = trim_coords_bed
        }
    }

    call read_utils.merge_and_reheader_bams as merge_align_to_ref {
        input:
            in_bams             = ivar_trim.aligned_trimmed_bam,
            sample_name         = sample_name,
            out_basename        = "${sample_name}.align_to_ref.trimmed"
    }

    call reports.plot_coverage as plot_ref_coverage {
        input:
            aligned_reads_bam   = merge_align_to_ref.out_bam,
            sample_name         = sample_name
    }

    call assembly.refine_assembly_with_aligned_reads as call_consensus {
        input:
            reference_fasta   = reference_fasta,
            reads_aligned_bam = merge_align_to_ref.out_bam,
            sample_name       = sample_name
    }

    scatter(reads_unmapped_bam in reads_unmapped_bams) {
        call assembly.align_reads as align_to_self {
            input:
                reference_fasta    = call_consensus.refined_assembly_fasta,
                reads_unmapped_bam = reads_unmapped_bam,
                novocraft_license  = novocraft_license,
                skip_mark_dupes    = skip_mark_dupes,
                aligner            = aligner,
                aligner_options    = align_to_self_options[aligner]
        }
    }

    call read_utils.merge_and_reheader_bams as merge_align_to_self {
        input:
            in_bams             = align_to_self.aligned_only_reads_bam,
            sample_name         = sample_name,
            out_basename        = "${sample_name}.merge_align_to_self"
    }

    call reports.plot_coverage as plot_self_coverage {
        input:
            aligned_reads_bam   = merge_align_to_self.out_bam,
            sample_name         = sample_name
    }

    output {
        File   assembly_fasta               = call_consensus.refined_assembly_fasta
        File   align_to_ref_variants_vcf_gz = call_consensus.sites_vcf_gz
        Int    assembly_length              = call_consensus.assembly_length
        Int    assembly_length_unambiguous  = call_consensus.assembly_length_unambiguous
        Int    reference_genome_length      = plot_ref_coverage.assembly_length
        Float  assembly_mean_coverage       = plot_ref_coverage.mean_coverage

        Array[File]   align_to_ref_per_input_aligned_flagstat = align_to_ref.aligned_bam_flagstat
        Array[Int]    align_to_ref_per_input_reads_provided   = align_to_ref.reads_provided
        Array[Int]    align_to_ref_per_input_reads_aligned    = align_to_ref.reads_aligned
        Array[File]   align_to_ref_per_input_fastqc = align_to_ref.aligned_only_reads_fastqc

        File   align_to_ref_merged_aligned_trimmed_only_bam = merge_align_to_ref.out_bam
        File   align_to_ref_merged_coverage_plot            = plot_ref_coverage.coverage_plot
        File   align_to_ref_merged_coverage_tsv             = plot_ref_coverage.coverage_tsv
        Int    align_to_ref_merged_reads_aligned            = plot_ref_coverage.reads_aligned
        Int    align_to_ref_merged_read_pairs_aligned       = plot_ref_coverage.read_pairs_aligned
        Float  align_to_ref_merged_bases_aligned            = plot_ref_coverage.bases_aligned

        File   align_to_self_merged_aligned_only_bam   = merge_align_to_self.out_bam
        File   align_to_self_merged_coverage_plot      = plot_self_coverage.coverage_plot
        File   align_to_self_merged_coverage_tsv       = plot_self_coverage.coverage_tsv
        Int    align_to_self_merged_reads_aligned      = plot_self_coverage.reads_aligned
        Int    align_to_self_merged_read_pairs_aligned = plot_self_coverage.read_pairs_aligned
        Float  align_to_self_merged_bases_aligned      = plot_self_coverage.bases_aligned
        Float  align_to_self_merged_mean_coverage            = plot_self_coverage.mean_coverage

        String align_to_ref_viral_core_version = align_to_ref.viralngs_version[0]
        String ivar_version                    = ivar_trim.ivar_version[0]
        String viral_assemble_version          = call_consensus.viralngs_version
    }

}
