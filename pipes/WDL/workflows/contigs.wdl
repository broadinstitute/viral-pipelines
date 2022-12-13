version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_taxon_filter.wdl" as taxon_filter
import "../tasks/tasks_assembly.wdl" as assembly

workflow contigs {
    input {
        File        reads_unmapped_bam

        Array[File] deplete_bmtaggerDbs = []
        Array[File] deplete_blastDbs = []
        Array[File] deplete_bwaDbs = []
    }

    if(length(deplete_bmtaggerDbs) + length(deplete_blastDbs) + length(deplete_bwaDbs) > 0) {
        call taxon_filter.deplete_taxa as deplete {
          input:
            raw_reads_unmapped_bam = reads_unmapped_bam,
            bmtaggerDbs            = deplete_bmtaggerDbs,
            blastDbs               = deplete_blastDbs,
            bwaDbs                 = deplete_bwaDbs
        }
    }

    call read_utils.rmdup_ubam {
        input:
            reads_unmapped_bam = select_first([deplete.cleaned_bam, reads_unmapped_bam])
    }

    call assembly.assemble as spades {
        input:
            reads_unmapped_bam = rmdup_ubam.dedup_bam
    }

    # TO DO: taxonomic classification of contigs

    output {
        File?   cleaned_reads_unaligned_bam     = deplete.cleaned_bam
        File    deduplicated_reads_unaligned    = rmdup_ubam.dedup_bam
        File    contigs_fastas                  = spades.contigs_fasta
        
        Int?    read_counts_raw                 = deplete.depletion_read_count_pre
        Int?    read_counts_depleted            = deplete.depletion_read_count_post
        Int     read_counts_dedup               = rmdup_ubam.dedup_read_count_post
        Int     read_counts_prespades_subsample = spades.subsample_read_count
        
        String? deplete_viral_classify_version  = deplete.viralngs_version
        String  spades_viral_assemble_version   = spades.viralngs_version
   }
}
