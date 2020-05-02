version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_assembly.wdl" as assembly

workflow contigs {

    call taxon_filter.deplete_taxa as deplete

    call assembly.assemble as spades {
        input:
            assembler = "spades",
            reads_unmapped_bam = deplete.cleaned_bam
    }

    # TO DO: taxonomic classification of contigs

    output {
        File   cleaned_reads_unaligned_bams    = deplete.cleaned_bam
        File   contigs_fastas                  = spades.contigs_fasta

        Int    read_counts_raw                 = deplete.depletion_read_count_pre
        Int    read_counts_depleted            = deplete.depletion_read_count_post
        Int    read_counts_prespades_subsample = spades.subsample_read_count

        String deplete_viral_classify_version  = deplete.viralngs_version
        String spades_viral_assemble_version   = spades.viralngs_version
   }
}
