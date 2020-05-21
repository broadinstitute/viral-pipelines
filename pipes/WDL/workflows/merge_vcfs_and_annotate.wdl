version 1.0

import "../tasks/tasks_interhost.wdl" as interhost
import "../tasks/tasks_intrahost.wdl" as intrahost

workflow merge_vcfs_and_annotate {
    meta {
        description: "Merge VCFs emitted by GATK UnifiedGenotyper and annotate with snpEff."
    }

    input {
        File reference_fasta
    }

    parameter_meta {
        reference_fasta: {
          description: "Reference genome, all segments/chromosomes in one fasta file. Headers must be Genbank accessions.",
          patterns: ["*.fasta"]
        }
    }

    call interhost.merge_vcfs_gatk as merge_vcfs {
        input:
            ref_fasta = reference_fasta
    }
    call intrahost.annotate_vcf_snpeff as annotate_vcf {
        input:
            ref_fasta = reference_fasta,
            in_vcf    = merge_vcfs.merged_vcf_gz
    }
    output {
        File merged_vcf_gz           = merge_vcfs.merged_vcf_gz
        File merged_vcf_gz_tbi       = merge_vcfs.merged_vcf_gz_tbi
        File merged_annot_vcf_gz     = annotate_vcf.annot_vcf_gz
        File merged_annot_vcf_gz_tbi = annotate_vcf.annot_vcf_gz_tbi
        File merged_annot_txt_gz     = annotate_vcf.annot_txt_gz
    }
}
