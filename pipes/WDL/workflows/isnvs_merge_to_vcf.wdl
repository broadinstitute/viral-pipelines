version 1.0

import "../tasks/tasks_interhost.wdl" as interhost
import "../tasks/tasks_intrahost.wdl" as tasks_intrahost

workflow isnvs_merge_to_vcf {
    input {
        File          reference_fasta
        Array[File]+  assemblies_fasta     # one per genome
    }

    call interhost.multi_align_mafft_ref as mafft {
        input:
            reference_fasta  = reference_fasta,
            assemblies_fasta = assemblies_fasta
    }

    call tasks_intrahost.isnvs_vcf {
        input:
            perSegmentMultiAlignments = mafft.alignments_by_chr,
            reference_fasta           = reference_fasta
    }

    output {
        Array[File] alignments_by_chr   = mafft.alignments_by_chr
        File        isnvs_plain_vcf     = isnvs_vcf.isnvs_vcf
        File        isnvs_plain_vcf_idx = isnvs_vcf.isnvs_vcf_idx
        File        isnvs_annot_vcf     = isnvs_vcf.isnvs_annot_vcf
        File        isnvs_annot_vcf_idx = isnvs_vcf.isnvs_annot_vcf_idx
        File        isnvs_annot_txt     = isnvs_vcf.isnvs_annot_txt
        String      viral_phylo_version = mafft.viralngs_version
    }
}
