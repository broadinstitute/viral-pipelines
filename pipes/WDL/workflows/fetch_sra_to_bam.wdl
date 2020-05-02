version 1.0

import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools

workflow multi_Fetch_SRA_to_BAM {
    input {
        Array[String]  SRR_accessions
    }
    scatter(SRA_ID in SRR_accessions) {
        call ncbi_tools.Fetch_SRA_to_BAM {
            input:
                SRA_ID = SRA_ID
        }
    }
    output {
        Array[File]   reads_ubam = Fetch_SRA_to_BAM.reads_ubam
        Array[String] sequencing_center = Fetch_SRA_to_BAM.sequencing_center
        Array[String] sequencing_platform = Fetch_SRA_to_BAM.sequencing_platform
        Array[String] sequencing_platform_model = Fetch_SRA_to_BAM.sequencing_platform_model
        Array[String] biosample_accession = Fetch_SRA_to_BAM.biosample_accession
        Array[String] library_id = Fetch_SRA_to_BAM.library_id
        Array[String] run_date = Fetch_SRA_to_BAM.run_date
        Array[String] sample_collection_date = Fetch_SRA_to_BAM.sample_collection_date
        Array[String] sample_collected_by = Fetch_SRA_to_BAM.sample_collected_by
        Array[String] sample_strain = Fetch_SRA_to_BAM.sample_strain
        Array[String] sample_geo_loc = Fetch_SRA_to_BAM.sample_geo_loc
        Array[File]   sra_metadata = Fetch_SRA_to_BAM.sra_metadata
    }
}
