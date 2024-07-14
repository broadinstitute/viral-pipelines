


workflow fetch_sra_to_bam_bam_only {
    meta {
        description: "Retrieve reads from the NCBI Short Read Archive in unaligned BAM format."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    call terra.check_terra_env

    #if(check_terra_env.is_running_on_terra) {
    call ncbi_tools.Fetch_SRA_to_BAM_BAM_only {
        input:
            email_address = select_first([check_terra_env.user_email, ""])
    }
    #}
    #if(!check_terra_env.is_running_on_terra) {
    #    call ncbi_tools.Fetch_SRA_to_BAM
    #}

    output {
        File   reads_ubam                = Fetch_SRA_to_BAM_BAM_only.reads_ubam
    }
}
