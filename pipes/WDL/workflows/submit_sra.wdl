version 1.0

import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_terra.wdl" as terra

workflow submit_sra {
    meta {
        description: "Submit reads to SRA"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[File]  reads_bams
        String       flowcell_id

        File         ncbi_ftp_config_js
        File         sra_meta_tsv
        String       sra_bioproject
        String       sra_data_bucket_uri

        String       prod_test = "Production" # Production or Test
    }

    String prefix = "/~{prod_test}/~{flowcell_id}"


    call terra.gcs_copy as gcs_sra_dump_reads {
        input:
            infiles        = reads_bams,
            gcs_uri_prefix = "~{sra_data_bucket_uri}/~{flowcell_id}/"
    }

    call ncbi_tools.sra_tsv_to_xml {
        input:
            meta_submit_tsv  = sra_meta_tsv,
            config_js        = ncbi_ftp_config_js,
            bioproject       = sra_bioproject,
            data_bucket_uri  = "~{sra_data_bucket_uri}/~{flowcell_id}"
    }
    call ncbi_tools.ncbi_sftp_upload as sra_upload {
        input:
            config_js        = ncbi_ftp_config_js,
            submission_xml   = sra_tsv_to_xml.submission_xml,
            additional_files = [],
            target_path      = "~{prefix}/sra",
            wait_for         = "1"
    }

    output {
        File          sra_xml            = sra_tsv_to_xml.submission_xml
        Array[File]   sra_response       = sra_upload.reports_xmls
    }
}