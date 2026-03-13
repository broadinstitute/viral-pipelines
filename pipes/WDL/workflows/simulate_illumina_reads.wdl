version 1.1

import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_assembly.wdl" as assembly

workflow simulate_illumina_reads {

    meta {
        description: "Generate synthetic Illumina read sets for testing using wgsim. Takes a space-separated string of colon-separated pairs where each pair consists of a GenBank accession and a coverage value (e.g., 'KJ660346.2:12.5x NC_004296.1:0.9X'), downloads the sequences, and simulates Illumina reads."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    parameter_meta {
        accession_coverage_string: {
            description: "Space-separated string of colon-separated pairs where each pair consists of a GenBank accession and a coverage value (e.g., 'KJ660346.2:12.5x NC_004296.1:0.9X').",
            category: "required"
        }
        out_basename: {
            description: "Base name for output files.",
            category: "common"
        }
    }

    input {
        String accession_coverage_string
        String out_basename = "simulated_reads"
    }

    call ncbi.download_fasta_from_accession_string {
        input:
            accession_string = accession_coverage_string,
            out_prefix       = out_basename
    }

    call assembly.wgsim as simulate_reads {
        input:
            coverage_string = accession_coverage_string,
            reference_fasta = download_fasta_from_accession_string.sequences_fasta,
            out_basename    = out_basename
    }

    output {
        File   simulated_reads_bam = simulate_reads.simulated_reads_bam
        Int    read_count          = simulate_reads.read_count
        File   reference_sequences = download_fasta_from_accession_string.sequences_fasta
        String viralngs_version    = simulate_reads.viralngs_version
    }
}
