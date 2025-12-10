version 1.0

import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow pairwise_distances {
    meta {
        description: "Download genomes by accession, compute pairwise distances with skani, and optionally generate MSA and tree"
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    parameter_meta {
        accessions_file: {
            description: "Text file with GenBank accession numbers for single segment genomes, one per line",
            category: "required"
        }
        emailAddress: {
            description: "Email address for NCBI Entrez identification",
            category: "required"
        }
        make_tree: {
            description: "If true, generate multiple sequence alignment and phylogenetic tree in addition to pairwise distances",
            category: "common"
        }
        out_prefix: {
            description: "Prefix for output file names",
            category: "common"
        }
    }

    input {
        File    accessions_file
        String  emailAddress
        Boolean make_tree = true
        String  out_prefix = "pairwise_dist"
    }

    # Read accessions from file (one per line)
    Array[String] accessions = read_lines(accessions_file)

    # Download sequences from GenBank
    call ncbi.download_fasta {
        input:
            out_prefix   = out_prefix,
            accessions   = accessions,
            emailAddress = emailAddress
    }

    # Compute pairwise distances with skani triangle
    call assembly.skani_triangle {
        input:
            sequences_fasta = download_fasta.sequences_fasta,
            out_basename    = out_prefix
    }

    # Optional: create MSA and tree
    # Note: Task inputs like ref_fasta, large, memsavetree for mafft_one_chr
    # are left unbound and will be exposed as optional nested inputs
    if (make_tree) {
        call nextstrain.mafft_one_chr {
            input:
                sequences = download_fasta.sequences_fasta,
                basename  = out_prefix + "_msa"
                # ref_fasta is intentionally unbound - user can provide via nested inputs
                # large and memsavetree can be set via nested inputs for 15k genomes
        }

        call nextstrain.draft_augur_tree {
            input:
                msa_or_vcf = mafft_one_chr.aligned_sequences
                # method, substitution_model, cpus, machine_mem_gb intentionally unbound
        }
    }

    output {
        File  sequences_fasta     = download_fasta.sequences_fasta
        File  pairwise_dist_skani = skani_triangle.skani_triangle_tsv
        File? msa_fasta           = mafft_one_chr.aligned_sequences
        File? iqtree_nwk          = draft_augur_tree.aligned_tree
    }
}
