version 1.0

import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_utils.wdl" as utils

workflow submit_biosample {
    meta {
        description: "Register samples with NCBI BioSample. Return attributes table."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File   biosample_submit_tsv
        File   ftp_config_js
        String prod_test = "Production" # Production or Test
    }

    call utils.md5sum {
        input:
            in_file = biosample_submit_tsv
    }

    # see if anything already exists in NCBI
    call ncbi_tools.biosample_tsv_filter_preexisting {
        input:
            meta_submit_tsv = biosample_submit_tsv,
            out_basename = basename(biosample_submit_tsv, '.tsv')
    }

    # register anything that isn't already in NCBI
    if (biosample_tsv_filter_preexisting.num_not_found > 0) {
        call ncbi_tools.biosample_submit_tsv_ftp_upload {
            input:
                meta_submit_tsv = biosample_tsv_filter_preexisting.meta_unsubmitted_tsv,
                config_js = ftp_config_js,
                target_path = "/~{prod_test}/biosample/~{basename(biosample_submit_tsv, '.tsv')}/~{md5sum.md5}"
        }
    }

    # merge all results and attributes
    call utils.tsv_join {
        input:
            input_tsvs = select_all([
                biosample_tsv_filter_preexisting.biosample_attributes_tsv,
                biosample_submit_tsv_ftp_upload.attributes_tsv,
                biosample_submit_tsv
            ]),
            id_col = "isolate",
            out_basename = basename(biosample_submit_tsv, '.tsv') + "-attributes"
    }

    output {
        File  biosample_attributes = tsv_join.out_tsv
    }
}
