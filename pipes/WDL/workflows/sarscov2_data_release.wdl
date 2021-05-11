version 1.0

import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_utils.wdl" as utils

workflow sarscov2_data_release {
    meta {
        description: "Submit data bundles to databases and repositories"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        String       flowcell_id

        File         ncbi_ftp_config_js
        File         genbank_xml
        File         genbank_zip
        #File         sra_meta_tsv

        #File         gisaid_csv
        #File         gisaid_fasta

        File?        cdc_s3_credentials
        File?        cdc_passing_fasta
        File?        cdc_final_metadata
        File?        cdc_cumulative_metadata
        Array[File]  cdc_aligned_trimmed_bams
        String?      cdc_s3_uri

        String       ftp_path_prefix = basename(genbank_zip, ".zip")
        String       prod_test = "Production" # Production or Test
    }

    String prefix = "/~{prod_test}/~{ftp_path_prefix}"

    call utils.today {
        input: timezone = "America/New_York"  # CDC is based in Atlanta
    }

    # publish to NCBI Genbank
    call ncbi_tools.ncbi_ftp_upload as genbank {
        input:
            config_js      = ncbi_ftp_config_js,
            submit_files   = [genbank_xml, genbank_zip],
            target_path    = "~{prefix}/genbank",
            wait_for       = "1"
    }

    # to do: Asymmetrik to impement an SRA tsv->xml conversion

    # to do: viral to implement a GISAID CLI upload task

    # deliver to CDC
    if (defined(cdc_s3_credentials)) {
        String s3_prefix = "~{cdc_s3_uri}/~{today.date}/~{flowcell_id}"
        call utils.make_empty_file as upload_complete {
            input:
                out_filename = "uploadcomplete.txt"
        }
        if (defined(cdc_cumulative_metadata)) {
            call utils.rename_file as cumulative_meta_tsv {
                input:
                    infile = select_first([cdc_cumulative_metadata]),
                    out_filename = "metadata-cumulative-~{today.date}.txt"
            }
            call utils.s3_copy as s3_cdc_dump_cumulative {
                input:
                    infiles         = [cumulative_meta_tsv.out],
                    s3_uri_prefix   = "~{cdc_s3_uri}/",
                    aws_credentials = select_first([cdc_s3_credentials])
            }
        }
        call utils.s3_copy as s3_cdc_dump_meta {
            input:
                infiles         = select_all([cdc_final_metadata, cdc_passing_fasta, upload_complete.out]),
                s3_uri_prefix   = "~{s3_prefix}/",
                aws_credentials = select_first([cdc_s3_credentials])
        }
        call utils.s3_copy as s3_cdc_dump_reads {
            input:
                infiles         = cdc_aligned_trimmed_bams,
                s3_uri_prefix   = "~{s3_prefix}/rawfiles/",
                aws_credentials = select_first([cdc_s3_credentials])
        }
    }

    output {
        Array[File]    genbank_response   = genbank.reports_xmls
    }
}
