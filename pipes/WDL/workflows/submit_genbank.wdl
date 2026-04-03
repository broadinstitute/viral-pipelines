version 1.1

import "../tasks/tasks_utils.wdl" as utils
import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools

workflow submit_genbank {
    meta {
        description: "Submit FTP-eligible genomes to NCBI Genbank (currently only flu A/B/C and SARS-CoV-2)"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]  submission_files # paired xml and zip files with same basenames
        String       batch_id

        File         ncbi_ftp_config_js
        String       prod_test = "Production" # Production or Test
    }

    String prefix = "/~{prod_test}/~{batch_id}"

    call utils.pair_files_by_basename {
        input:
            left_ext  = "xml",
            right_ext = "zip",
            files     = submission_files
    }

    scatter(file_pair in pair_files_by_basename.file_pairs) {
        call ncbi_tools.ncbi_sftp_upload as genbank_upload {
            input:
                config_js        = ncbi_ftp_config_js,
                submission_xml   = file_pair.left,
                additional_files = [file_pair.right],
                target_path      = "~{prefix}/genbank/~{basename(file_pair.right)}",
                wait_for         = "1"
        }
    }

    output {
        Array[File]    genbank_response   = flatten(genbank_upload.reports_xmls)
    }
}