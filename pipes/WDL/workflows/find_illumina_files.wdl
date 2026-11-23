version 1.1

#DX_SKIP_WORKFLOW

import "../tasks/tasks_terra.wdl" as terra

workflow find_illumina_files {
    meta {
        description: "Utility workflow to discover and structure Illumina run files from a GCS bucket directory. This workflow uses gcloud storage commands and is only compatible with Google Cloud Platform."
        author: "Broad Viral Genomics"
        email: "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    call terra.find_illumina_files_in_directory

    output {
        String               runinfo_xml           = find_illumina_files_in_directory.runinfo_xml
        Array[String]        fastqs                = find_illumina_files_in_directory.fastqs
    }
}
