version 1.0

import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_sarscov2.wdl" as sarscov2


workflow sarscov2_biosample_load {
    meta {
        description: "Load Broad CRSP metadata and register samples with NCBI BioSample. Return attributes table, id map, etc."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        File?  sample_meta_crsp
        File   id_salt
        File?  biosample_submit_tsv
        String bioproject
        File   ftp_config_js
        String ftp_target_path
    }

    if(!defined(biosample_submit_tsv)) {
        call sarscov2.crsp_meta_etl {
            input:
                sample_meta_crsp = select_first([sample_meta_crsp]),
                salt = read_string(id_salt),
                bioproject = bioproject
        }
    }

    call ncbi_tools.biosample_submit_tsv_ftp_upload {
        input:
            meta_submit_tsv = select_first([biosample_submit_tsv, crsp_meta_etl.biosample_submit_tsv]),
            config_js = ftp_config_js,
            target_path = ftp_target_path
    } 

    output {
        File           biosample_attributes = biosample_submit_tsv_ftp_upload.attributes_tsv
        File?          id_map_tsv           = crsp_meta_etl.collab_ids_tsv
        File?          collab_ids_tsv       = crsp_meta_etl.collab_ids_tsv
        Array[String]  collab_ids_addcols   = select_first([crsp_meta_etl.collab_ids_addcols, []])
    }

}
