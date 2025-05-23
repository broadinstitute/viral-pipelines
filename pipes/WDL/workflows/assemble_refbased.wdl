version 1.0

import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_intrahost.wdl" as intrahost
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_read_utils.wdl" as read_utils

workflow assemble_refbased {

    meta {
        description: "Reference-based microbial consensus calling. Aligns NGS reads to a singular reference genome, calls a new consensus sequence, and emits: new assembly, reads aligned to provided reference, reads aligned to new assembly, various figures of merit, plots, and QC metrics. The user may provide unaligned reads spread across multiple input files and this workflow will parallelize alignment per input file before merging results prior to consensus calling."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
        # https://wdl-aid.readthedocs.io/en/latest/usage.html -- can exclude certain inputs from documentation; attempt to do this on with docker nominated variables 
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
        skip_mark_dupes: {
            description: "skip Picard MarkDuplicates step after alignment. This is recommended to be set to true for PCR amplicon based data. (Default: false)"
        }
        trim_coords_bed: {
            description: "optional primers to trim in reference coordinate space (0-based BED format)",
            patterns: ["*.bed"],
            category: "common"
        }
        min_coverage: {
            description: "Minimum read coverage required to call a position unambiguous."
        }
        major_cutoff: {
            description: "If the major allele is present at a frequency higher than this cutoff, we will call an unambiguous base at that position.  If it is equal to or below this cutoff, we will call an ambiguous base representing all possible alleles at that position."
        }
        assembly_fasta: { 
            description: "The new assembly / consensus sequence for this sample" 
        }
        align_to_ref_variants_vcf_gz: { 
            description: "All variants in the input reads against the original reference genome. This VCF file is used to create the assembly_fasta" 
        }
        assembly_length: { 
            description: "The length of the sequence described in assembly_fasta, inclusive of any uncovered regions denoted by Ns" 
        }
        assembly_length_unambiguous: { 
            description: "The number of called consensus bases in assembly_fasta (excludes regions of the genome that lack read coverage)" 
        }
    }

    input {
        Array[File]+ reads_unmapped_bams
        File         reference_fasta
        String       sample_name = basename(reads_unmapped_bams[0], '.bam')
        String?      sample_original_name

        String       aligner="minimap2"
        File?        novocraft_license
        Int          min_coverage=3
        Float        major_cutoff=0.75
        Boolean      skip_mark_dupes=false
        File?        trim_coords_bed

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
                aligned_bam     = align_to_ref.aligned_only_reads_bam,
                trim_coords_bed = trim_coords_bed
        }
        Map[String,String] ivar_stats = {
            'file': basename(reads_unmapped_bam, '.bam'),
            'trim_percent': ivar_trim.primer_trimmed_read_percent,
            'trim_count':   ivar_trim.primer_trimmed_read_count
        }
        Array[String] ivar_stats_row = [ivar_stats['file'], ivar_stats['trim_percent'], ivar_stats['trim_count']]
    }

    if(length(reads_unmapped_bams)>1) {
        call read_utils.merge_and_reheader_bams as merge_align_to_ref {
            input:
                in_bams      = ivar_trim.aligned_trimmed_bam,
                sample_name  = sample_name,
                out_basename = "~{sample_name}.align_to_ref.trimmed"
        }
    }
    File aligned_trimmed_bam = select_first([merge_align_to_ref.out_bam, ivar_trim.aligned_trimmed_bam[0]])

    call intrahost.lofreq as isnvs_ref {
        input:
            reference_fasta = reference_fasta,
            aligned_bam     = aligned_trimmed_bam
    }

    call reports.alignment_metrics {
        input:
            aligned_bam = aligned_trimmed_bam,
            ref_fasta   = reference_fasta,
            primers_bed = trim_coords_bed,
            min_coverage = min_coverage
    }

    call assembly.run_discordance {
        input:
            reads_aligned_bam = aligned_trimmed_bam,
            reference_fasta   = reference_fasta,
            min_coverage      = min_coverage+1,
            out_basename      = sample_name
    }

    call reports.plot_coverage as plot_ref_coverage {
        input:
            aligned_reads_bam = aligned_trimmed_bam,
            sample_name       = sample_name
    }

    call assembly.refine_assembly_with_aligned_reads as call_consensus {
        input:
            reference_fasta   = reference_fasta,
            reads_aligned_bam = aligned_trimmed_bam,
            min_coverage      = min_coverage,
            major_cutoff      = major_cutoff,
            out_basename      = sample_name,
            sample_name       = select_first([sample_original_name, sample_name])
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

    if(length(reads_unmapped_bams)>1) {
        call read_utils.merge_and_reheader_bams as merge_align_to_self {
            input:
                in_bams      = align_to_self.aligned_only_reads_bam,
                sample_name  = sample_name,
                out_basename = "~{sample_name}.merge_align_to_self"
        }
    }
    File aligned_self_bam = select_first([merge_align_to_self.out_bam, align_to_self.aligned_only_reads_bam[0]])

    call intrahost.lofreq as isnvs_self {
        input:
            reference_fasta = call_consensus.refined_assembly_fasta,
            aligned_bam     = aligned_self_bam
    }

    call reports.plot_coverage as plot_self_coverage {
        input:
            aligned_reads_bam = aligned_self_bam,
            sample_name       = sample_name
    }

    output {
        File        assembly_fasta                               = call_consensus.refined_assembly_fasta
        File        align_to_ref_variants_vcf_gz                 = call_consensus.sites_vcf_gz
        Int         assembly_length                              = call_consensus.assembly_length
        Int         assembly_length_unambiguous                  = call_consensus.assembly_length_unambiguous
        Int         reference_genome_length                      = align_to_ref.reference_length[0]
        Float       assembly_mean_coverage                       = plot_ref_coverage.mean_coverage
        
        Int         dist_to_ref_snps                             = call_consensus.dist_to_ref_snps
        Int         dist_to_ref_indels                           = call_consensus.dist_to_ref_indels
        
        Array[Int]                primer_trimmed_read_count      = ivar_trim.primer_trimmed_read_count
        Array[Float]              primer_trimmed_read_percent    = ivar_trim.primer_trimmed_read_percent
        Array[Map[String,String]] ivar_trim_stats                = ivar_stats
        Array[Array[String]]      ivar_trim_stats_tsv            = ivar_stats_row
        
        Int         replicate_concordant_sites                   = run_discordance.concordant_sites
        Int         replicate_discordant_snps                    = run_discordance.discordant_snps
        Int         replicate_discordant_indels                  = run_discordance.discordant_indels
        Int         num_read_groups                              = run_discordance.num_read_groups
        Int         num_libraries                                = run_discordance.num_libraries
        File        replicate_discordant_vcf                     = run_discordance.discordant_sites_vcf
        
        Array[File] align_to_ref_per_input_aligned_flagstat      = align_to_ref.aligned_bam_flagstat
        Array[Int]  align_to_ref_per_input_reads_provided        = align_to_ref.reads_provided
        Array[Int]  align_to_ref_per_input_reads_aligned         = align_to_ref.reads_aligned

        File        align_to_ref_merged_aligned_trimmed_only_bam = aligned_trimmed_bam
        File        align_to_ref_fastqc                          = select_first([merge_align_to_ref.fastqc, align_to_ref.aligned_only_reads_fastqc[0]])
        File        align_to_ref_merged_coverage_plot            = plot_ref_coverage.coverage_plot
        File        align_to_ref_merged_coverage_tsv             = plot_ref_coverage.coverage_tsv
        Int         align_to_ref_merged_reads_aligned            = plot_ref_coverage.reads_aligned
        Int         align_to_ref_merged_read_pairs_aligned       = plot_ref_coverage.read_pairs_aligned
        Int         align_to_ref_merged_bases_aligned            = floor(plot_ref_coverage.bases_aligned)
        File        align_to_ref_isnvs_vcf                       = isnvs_ref.report_vcf
        
        File        picard_metrics_wgs                           = alignment_metrics.wgs_metrics
        File        picard_metrics_alignment                     = alignment_metrics.alignment_metrics
        File        picard_metrics_insert_size                   = alignment_metrics.insert_size_metrics
        File        samtools_ampliconstats                       = alignment_metrics.amplicon_stats
        File        samtools_ampliconstats_parsed                = alignment_metrics.amplicon_stats_parsed

        Array[File] align_to_self_merged_aligned_and_unaligned_bam = align_to_self.aligned_bam

        File        align_to_self_merged_aligned_only_bam        = aligned_self_bam
        File        align_to_self_merged_coverage_plot           = plot_self_coverage.coverage_plot
        File        align_to_self_merged_coverage_tsv            = plot_self_coverage.coverage_tsv
        Int         align_to_self_merged_reads_aligned           = plot_self_coverage.reads_aligned
        Int         align_to_self_merged_read_pairs_aligned      = plot_self_coverage.read_pairs_aligned
        Int         align_to_self_merged_bases_aligned           = floor(plot_self_coverage.bases_aligned)
        Float       align_to_self_merged_mean_coverage           = plot_self_coverage.mean_coverage
        File        align_to_self_isnvs_vcf                      = isnvs_self.report_vcf
        
        String      assembly_method = "viral-ngs/assemble_refbased"
        String      align_to_ref_viral_core_version              = align_to_ref.viralngs_version[0]
        String      ivar_version                                 = ivar_trim.ivar_version[0]
        String      viral_assemble_version                       = call_consensus.viralngs_version
    }

}
