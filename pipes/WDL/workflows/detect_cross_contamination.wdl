version 1.1

import "../tasks/tasks_intrahost.wdl" as intrahost
import "../tasks/tasks_read_utils.wdl" as read_utils

workflow detect_cross_contamination {
    meta {
        description: "Detect cross-contamination between samples using consensus-level and sub-consensus variation."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    parameter_meta {
        aligned_bams: {
            description: "reads aligned to the sequence within reference_fasta.",
            patterns: ["*.bam"]
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
        Array[File]+ aligned_bams
        Array[File]+ genome_fastas
        File         reference_fasta
    }

    scatter(sample in zip(aligned_bams,genome_fastas)) {
        call intrahost.lofreq as lofreq {
            input:
                aligned_bam = sample.left,
                reference_fasta = reference_fasta
        }
        
        call read_utils.read_depths as depth {
            input:
                aligned_bam = sample.left
        }
        
        # pair the lofreq and read depth outputs with the corresponding (input) consensus
        # fasta since we do not actually use the consensus fasta in the scattered task
        # (and thus cannot access it by task_name.input_name reference outside the scatter)
        # and since the WDL 1.0 spec is not explicit about whether order is preserved
        # during the gather phase (executor implementations may vary)
        # https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md
        Array[File] vcf_genome_read_depth_triplets = [
            lofreq.report_vcf,
            sample.right,
            depth.read_depths
        ]
    }

    Array[Array[File]] vcfs_and_genomes = transpose(vcf_genome_read_depth_triplets)
    call intrahost.polyphonia_detect_cross_contamination as detect_cross_contam {
        # take scatter-gathered array of [(vcf1,fasta1),(vcf2,fasta2),(vcf3,fasta3)]
        # and transpose to [[vcf1,vcf2,vcf3],[fasta1,fasta2,fasta3]]
        input:
            lofreq_vcfs     = vcfs_and_genomes[0], # vcfs
            genome_fastas   = vcfs_and_genomes[1], # fastas
            read_depths     = vcfs_and_genomes[2], # read depth tables
            reference_fasta = reference_fasta
    }

    output {
        File        contamination_report  = detect_cross_contam.report
        Array[File] lofreq_vcfs           = lofreq.report_vcf
        Array[File] read_depths           = depth.read_depths
        Array[File] contamination_figures = detect_cross_contam.figures
        String      lofreq_version        = lofreq.lofreq_version[0]
        # commented out until polyphonia can report its own version
        #String      polyphonia_version    = detect_cross_contam.polyphonia_version
    }
}
