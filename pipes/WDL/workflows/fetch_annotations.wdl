version 1.0

import "../tasks/tasks_ncbi.wdl" as ncbi

workflow fetch_annotations {
    call ncbi.download_annotations
    output {
        File        combined_fasta      = download_annotations.combined_fasta
        Array[File] genomes_fasta       = download_annotations.genomes_fasta
        Array[File] features_tbl        = download_annotations.features_tbl
        String      viral_phylo_version = download_annotations.viralngs_version
    }
}
