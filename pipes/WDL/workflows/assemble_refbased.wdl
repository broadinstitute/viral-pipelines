import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_read_utils.wdl" as read_utils

workflow assemble_refbased {

    String          sample_name
    Array[File]+    reads_unmapped_bams
    File            reference_fasta

    File?           novocraft_license
    Boolean?        skip_mark_dupes=false

    scatter(reads_unmapped_bam in reads_unmapped_bams) {
        call assembly.align_reads as align_to_ref {
            input:
                reference_fasta    = reference_fasta,
                reads_unmapped_bam = reads_unmapped_bam,
                novocraft_license  = novocraft_license,
                skip_mark_dupes    = skip_mark_dupes,
                aligner_options    = "-r Random -l 40 -g 40 -x 20 -t 501 -k"
                ## (for bwa) -- aligner_options = "-k 12 -B 1"
                ## (for novoalign) -- aligner_options = "-r Random -l 40 -g 40 -x 20 -t 501 -k"
        }
    }

    call read_utils.merge_bams_one_sample as merge_align_to_ref {
        input:
            in_bams             = align_to_ref.aligned_only_reads_bam,
            optional_reads_bais = align_to_ref.aligned_only_reads_bam_idx,
            sample_name         = sample_name,
            out_basename        = "${sample_name}.align_to_ref"
    }

    call reports.plot_coverage as plot_ref_coverage {
        input:
            aligned_reads_bam   = merge_align_to_ref.out_bam
    }

    call assembly.ivar_trim {
        input:
            aligned_bam = merge_align_to_ref.out_bam
    }

    call assembly.refine_assembly_with_aligned_reads as call_consensus {
        input:
            reference_fasta   = reference_fasta,
            reads_aligned_bam = ivar_trim.aligned_trimmed_bam,
            sample_name       = sample_name,
            novocraft_license = novocraft_license
    }

    scatter(reads_unmapped_bam in reads_unmapped_bams) {
        call assembly.align_reads as align_to_self {
            input:
                reference_fasta    = call_consensus.refined_assembly_fasta,
                reads_unmapped_bam = reads_unmapped_bam,
                novocraft_license  = novocraft_license,
                skip_mark_dupes    = skip_mark_dupes,
                aligner_options    = "-r Random -l 40 -g 40 -x 20 -t 100"
                ## (for bwa) -- aligner_options = "-k 12 -B 1"
                ## (for novoalign) -- aligner_options = "-r Random -l 40 -g 40 -x 20 -t 501 -k"
        }
    }

    call read_utils.merge_bams_one_sample as merge_align_to_self {
        input:
            in_bams             = align_to_self.aligned_only_reads_bam,
            optional_reads_bais = align_to_self.aligned_only_reads_bam_idx,
            sample_name         = sample_name,
            out_basename        = "${sample_name}.merge_align_to_self"
    }

    call reports.plot_coverage as plot_self_coverage {
        input:
            aligned_reads_bam   = merge_align_to_self.out_bam
    }
}
