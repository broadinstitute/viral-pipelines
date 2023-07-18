version 1.0

import "../tasks/tasks_ncbi.wdl" as ncbi

workflow fetch_fasta_for_genbank_accessions {
    meta {
        description: "Retrieve sequences from NCBI GenBank based on a list of accessions, and output a single fasta file containing the sequences."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    call ncbi.download_fasta

    output {
        File sequences_fasta = download_fasta.sequences_fasta
    }
}
