version 1.0

import "../tasks/tasks_interhost.wdl" as interhost
import "../tasks/tasks_ncbi.wdl" as ncbi

workflow genbank {

    input {
        File          reference_fasta
        Array[File]+  assemblies_fasta     # one per genome
    }

    call interhost.multi_align_mafft_ref as mafft {
        input:
            reference_fasta = reference_fasta,
            assemblies_fasta = assemblies_fasta
    }

    call ncbi.annot_transfer as annot {
        input:
            multi_aln_fasta = mafft.alignments_by_chr,
            reference_fasta = reference_fasta
    }
 
    call ncbi.prepare_genbank as prep_genbank {
        input:
            assemblies_fasta = assemblies_fasta,
            annotations_tbl = annot.transferred_feature_tables
    }

    output {
        Array[File] alignments_by_chr          = mafft.alignments_by_chr

        Array[File] transferred_feature_tables = annot.transferred_feature_tables

        Array[File] sequin_files               = prep_genbank.sequin_files
        Array[File] structured_comment_files   = prep_genbank.structured_comment_files
        Array[File] genbank_preview_files      = prep_genbank.genbank_preview_files
        Array[File] source_table_files         = prep_genbank.source_table_files
        Array[File] fasta_per_chr_files        = prep_genbank.fasta_per_chr_files
        Array[File] validation_files           = prep_genbank.validation_files
        File        errorSummary               = prep_genbank.errorSummary

        String      viral_phylo_version        = mafft.viralngs_version
    }

}
