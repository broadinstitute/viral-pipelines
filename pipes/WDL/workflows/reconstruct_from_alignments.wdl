version 1.1

import "../tasks/tasks_interhost.wdl" as interhost
import "../tasks/tasks_intrahost.wdl" as intrahost
import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_utils.wdl" as utils

workflow reconstruct_from_alignments {
    meta {
        description: "Infer disease transmission events from sequence (consensus + intrahost variation) data using the reconstructR tool"
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
    call utils.zcat {
        input:
            infiles     = assembly_fastas,
            output_name = "all_assemblies.fasta.gz"
    }
    call nextstrain.mafft_one_chr as mafft {
        input:
            sequences = zcat.combined,
            ref_fasta = ref_fasta,
            remove_reference = true,
            basename  = "mafft_msa.fasta"
    }

    # call iSNVs with lofreq and calculate coverage depths
    scatter(bam in aligned_trimmed_bams) {
        call read_utils.get_bam_samplename {
            input:
                bam = bam
        }
        call intrahost.lofreq as isnvs_ref {
            input:
                reference_fasta = ref_fasta,
                aligned_bam     = bam,
                out_basename    = get_bam_samplename.sample_name
        }
        call reports.plot_coverage as plot_ref_coverage {
            input:
                aligned_reads_bam = bam,
                sample_name       = get_bam_samplename.sample_name
        }
    }
    call reports.merge_coverage_per_position {
        input:
            coverage_tsvs = plot_ref_coverage.coverage_tsv,
            ref_fasta = ref_fasta
    }

    # call reconstructR
    call interhost.reconstructr {
        input:
            ref_fasta = ref_fasta,
            lofreq_vcfs = isnvs_ref.report_vcf,
            msa_fasta = mafft.aligned_sequences,
            depth_csv = merge_coverage_per_position.coverage_multi_sample_per_position_csv
    }

    output {
      File         msa_fasta                      = mafft.aligned_sequences
      Array[File]  lofreq_isnvs                   = isnvs_ref.report_vcf
      File         depth_csv                      = merge_coverage_per_position.coverage_multi_sample_per_position_csv
      File         reconstructr_tabulated_tsv_gz  = reconstructr.tabulated_tsv_gz
      File         reconstructr_deciphered_tsv_gz = reconstructr.deciphered_tsv_gz
    }
}