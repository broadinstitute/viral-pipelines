version 1.0

import "../tasks/tasks_interhost.wdl" as interhost

workflow merge_vcfs {
    meta {
        description: "Merge VCFs from multiple samples using GATK3."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }
    call interhost.merge_vcfs_gatk
    output {
        File merged_vcf_gz       = merge_vcfs_gatk.merged_vcf_gz
        File merged_vcf_gz_tbi   = merge_vcfs_gatk.merged_vcf_gz_tbi
    }
}
