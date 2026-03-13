version 1.1

import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_sarscov2.wdl" as sarscov2
import "../tasks/tasks_utils.wdl" as utils


workflow sarscov2_biosample_load {
    meta {
        description: "Load Broad CRSP metadata and register samples with NCBI BioSample. Return attributes table, id map, etc."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File?  sample_meta_crsp
        File   id_salt
        File?  biosample_submit_tsv
        String bioproject
        File   ftp_config_js
        String prod_test = "Production" # Production or Test
    }

    if(!defined(biosample_submit_tsv)) {
        call sarscov2.crsp_meta_etl {
            input:
                sample_meta_crsp = select_first([sample_meta_crsp]),
                salt = read_string(id_salt),
                bioproject = bioproject
        }
    }

    File meta_submit_tsv = select_first([biosample_submit_tsv, crsp_meta_etl.biosample_submit_tsv])

    call utils.md5sum {
        input:
            in_file = meta_submit_tsv
    }

    # see if anything already exists in NCBI
    call ncbi_tools.biosample_tsv_filter_preexisting {
        input:
            meta_submit_tsv = meta_submit_tsv,
            out_basename = basename(meta_submit_tsv, '.tsv')
    }

    # register anything that isn't already in NCBI
    if (biosample_tsv_filter_preexisting.num_not_found > 0) {
        call ncbi_tools.biosample_submit_tsv_ftp_upload {
            input:
                meta_submit_tsv = biosample_tsv_filter_preexisting.meta_unsubmitted_tsv,
                config_js = ftp_config_js,
                target_path = "/~{prod_test}/biosample/~{basename(meta_submit_tsv, '.tsv')}/~{md5sum.md5}"
        }
    }

    # merge all results and attributes
    call utils.tsv_join {
        input:
            input_tsvs = select_all([
                biosample_tsv_filter_preexisting.biosample_attributes_tsv,
                biosample_submit_tsv_ftp_upload.attributes_tsv,
                meta_submit_tsv
            ]),
            id_col = "isolate",
            out_basename = basename(meta_submit_tsv, '.tsv') + "-attributes"
    }

    output {
        File           biosample_attributes = tsv_join.out_tsv
        File?          id_map_tsv           = crsp_meta_etl.collab_ids_tsv
        File?          collab_ids_tsv       = crsp_meta_etl.collab_ids_tsv
        Array[String]  collab_ids_addcols   = select_first([crsp_meta_etl.collab_ids_addcols, []])
    }

}
