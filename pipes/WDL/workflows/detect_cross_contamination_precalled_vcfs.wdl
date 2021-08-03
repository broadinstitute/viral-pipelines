version 1.0

import "../tasks/tasks_intrahost.wdl" as intrahost

workflow detect_cross_contamination_precalled_vcfs {
    meta {
        description: "Detect cross-contamination between samples using consensus-level and sub-consensus variation, from consensus genomes and pre-called LoFreq vcf files."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    parameter_meta {
        lofreq_vcfs: {
            description: "per-sample variant call files (VCF) called by LoFreq against reference_fasta",
            patterns: ["*.vcf","*.vcf.gz"]
        }
        genome_fastas: {
            description: "consensus sequences, one per file in aligned_bams, in corresponding order",
            patterns: ["*.fasta"]
        }
        reference_fasta: {
            description: "Reference genome to which reads have be aligned, needed here for variant calling",
            patterns: ["*.fasta"]
        }

    }

    input {
        Array[File]+ lofreq_vcfs
        Array[File]+ genome_fastas
        File         reference_fasta
    }

    call intrahost.detect_cross_contamination as detect_cross_contam {
        input:
            lofreq_vcfs     = lofreq_vcfs,
            genome_fastas   = genome_fastas,
            reference_fasta = reference_fasta
    }

    output {
        File        contamination_report  = detect_cross_contam.report
        Array[File] contamination_figures = detect_cross_contam.figures
        # commented out until polyphonia can report its own version
        #String      polyphonia_version    = detect_cross_contam.polyphonia_version
    }
}
