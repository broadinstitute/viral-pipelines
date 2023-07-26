version 1.0

import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools

workflow fetch_sra_bams_for_genbank_accession {
    meta {
        description: "Retrieve reads from the NCBI Short Read Archive in unaligned BAM format, one file SRA run, associated with a given GenBank sequence accession."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    call ncbi_tools.fetch_sra_run_accessions_for_genbank_accession

    scatter(sra_accession in fetch_sra_run_accessions_for_genbank_accession.sra_accessions) {
        call ncbi_tools.Fetch_SRA_to_BAM {
            input:
                SRA_ID = sra_accession
        }
    }

    output {
        Array[File] sra_bams = Fetch_SRA_to_BAM.reads_ubam
    }
}
