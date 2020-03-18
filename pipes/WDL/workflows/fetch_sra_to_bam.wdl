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
}
