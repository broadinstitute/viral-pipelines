version 1.0

import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_utils.wdl" as utils

workflow genbank_gather {
    meta {
        description: "More here."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[String?]     genbank_file_manifest
        Array[File?]       genbank_submit_files
    }

    call ncbi.package_genbank_submissions {
        input:
            genbank_file_manifest   = select_all(genbank_file_manifest),
            genbank_submit_files    = select_all(genbank_submit_files)
    }

    output {
        Array[File] ftp_submission_files  = package_genbank_submissions.ftp_submission_files
        File? submit_sqns_clean_zip       = package_genbank_submissions.submit_sqns_clean_zip
        File? submit_sqns_warnings_zip    = package_genbank_submissions.submit_sqns_warnings_zip
        Int   num_sqns_clean              = package_genbank_submissions.num_sqns_clean
        Int   num_sqns_warnings           = package_genbank_submissions.num_sqns_warnings
        File? submit_fluA_clean_zip       = package_genbank_submissions.submit_fluA_clean_zip
        File? submit_fluA_warnings_zip    = package_genbank_submissions.submit_fluA_warnings_zip
        Int   num_fluA_clean              = package_genbank_submissions.num_fluA_clean
        Int   num_fluA_warnings           = package_genbank_submissions.num_fluA_warnings
        File? submit_fluB_clean_zip       = package_genbank_submissions.submit_fluB_clean_zip
        File? submit_fluB_warnings_zip    = package_genbank_submissions.submit_fluB_warnings_zip
        Int   num_fluB_clean              = package_genbank_submissions.num_fluB_clean
        Int   num_fluB_warnings           = package_genbank_submissions.num_fluB_warnings
        File? submit_fluC_clean_zip       = package_genbank_submissions.submit_fluC_clean_zip
        File? submit_fluC_warnings_zip    = package_genbank_submissions.submit_fluC_warnings_zip
        Int   num_fluC_clean              = package_genbank_submissions.num_fluC_clean
        Int   num_fluC_warnings           = package_genbank_submissions.num_fluC_warnings
        File? submit_sc2_clean_zip        = package_genbank_submissions.submit_sc2_clean_zip
        File? submit_sc2_warnings_zip     = package_genbank_submissions.submit_sc2_warnings_zip
        Int   num_sc2_clean               = package_genbank_submissions.num_sc2_clean
        Int   num_sc2_warnings            = package_genbank_submissions.num_sc2_warnings
        File? submit_noro_clean_zip       = package_genbank_submissions.submit_noro_clean_zip
        File? submit_noro_warnings_zip    = package_genbank_submissions.submit_noro_warnings_zip
        Int   num_noro_clean              = package_genbank_submissions.num_noro_clean
        Int   num_noro_warnings           = package_genbank_submissions.num_noro_warnings
        File? submit_dengue_clean_zip     = package_genbank_submissions.submit_dengue_clean_zip
        File? submit_dengue_warnings_zip  = package_genbank_submissions.submit_dengue_warnings_zip
        Int   num_dengue_clean            = package_genbank_submissions.num_dengue_clean
        Int   num_dengue_warnings         = package_genbank_submissions.num_dengue_warnings
    }
}