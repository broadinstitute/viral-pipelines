version 1.0

import "../tasks/tasks_interhost.wdl" as interhost
import "../tasks/tasks_intrahost.wdl" as intrahost
import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_utils.wdl" as utils

workflow reconstruct_from_alignments {
    meta {
        description: "TO DO"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]+   aligned_trimmed_bams
        Array[File]+   assembly_fastas
        File           ref_fasta
    }

    # create multiple sequence alignment of fastas
    call nextstrain.mafft_one_chr as mafft {
        input:
            sequences = assembly_fastas,
            ref_fasta = ref_fasta,
            basename  = "all_samples_aligned.fasta"
    }

    # call iSNVs with lofreq and calculate coverage depths
    scatter(bam in aligned_trimmed_bams) {
        call intrahost.lofreq as isnvs_ref {
            input:
                reference_fasta = ref_fasta,
                aligned_bam     = bam
        }
        call reports.plot_coverage as plot_ref_coverage {
            input:
                aligned_reads_bam = bam,
                sample_name       = basename(bam, '.bam')
        }
    }

    # TO DO: create depth_csv from plot_ref_coverage.coverage_tsv

    # call reconstructR
    call interhost.reconstructr {
        input:
            ref_fasta = ref_fasta,
            lofreq_vcfs = isnvs_ref.report_vcf,
            msa_fasta = mafft.aligned_sequences
        # TO DO: depth_csv = something computed above
    }

    output {
      File         msa_fasta                      = mafft.aligned_sequences
      Array[File]  lofreq_isnvs                   = isnvs_ref.report_vcf
      File         reconstructr_tabulated_tsv_gz  = reconstructr.tabulated_tsv_gz
      File         reconstructr_deciphered_tsv_gz = reconstructr.deciphered_tsv_gz
    }
}