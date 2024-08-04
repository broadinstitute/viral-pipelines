version 1.0

import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_terra.wdl" as terra

workflow fetch_sra_to_bam {
    meta {
        description: "Retrieve reads from the NCBI Short Read Archive in unaligned BAM format with relevant metadata encoded."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        String? email_address
    }

    call terra.check_terra_env

    String? user_email_address = select_first([email_address,check_terra_env.user_email, ""])

    #if(check_terra_env.is_running_on_terra) {
    call ncbi_tools.Fetch_SRA_to_BAM {
        input:
            email_address = user_email_address
    }
    #}
    #if(!check_terra_env.is_running_on_terra) {
    #    call ncbi_tools.Fetch_SRA_to_BAM
    #}

    output {
        File   reads_ubam                = Fetch_SRA_to_BAM.reads_ubam
        String sequencing_center         = Fetch_SRA_to_BAM.sequencing_center
        String sequencing_platform       = Fetch_SRA_to_BAM.sequencing_platform
        String sequencing_platform_model = Fetch_SRA_to_BAM.sequencing_platform_model
        String biosample_accession       = Fetch_SRA_to_BAM.biosample_accession
        String library_id                = Fetch_SRA_to_BAM.library_id
        String run_date                  = Fetch_SRA_to_BAM.run_date
        String sample_collection_date    = Fetch_SRA_to_BAM.sample_collection_date
        String sample_collected_by       = Fetch_SRA_to_BAM.sample_collected_by
        String sample_strain             = Fetch_SRA_to_BAM.sample_strain
        String sample_geo_loc            = Fetch_SRA_to_BAM.sample_geo_loc
        File   sra_metadata              = Fetch_SRA_to_BAM.sra_metadata
    }
}
