version 1.0

import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_terra.wdl" as terra

workflow submit_sra {
    meta {
        description: "Submit FTP-eligible genomes to NCBI Genbank (currently only flu A/B/C and SARS-CoV-2)"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        #Array[File]  submission_files
        String       batch_id

        File         genbank_xml
        File         genbank_zip

        File         ncbi_ftp_config_js
        String       prod_test = "Production" # Production or Test
    }

    String prefix = "/~{prod_test}/~{batch_id}"

    # TO DO: work in progress, not ready yet

    call ncbi_tools.ncbi_sftp_upload as genbank_upload {
        input:
            config_js        = ncbi_ftp_config_js,
            submission_xml   = genbank_xml,
            additional_files = [genbank_zip],
            target_path      = "~{prefix}/genbank",
            wait_for         = "1"
    }

    output {
        Array[File]    genbank_response   = genbank_upload.reports_xmls
    }
}