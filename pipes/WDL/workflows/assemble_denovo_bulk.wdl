import "tasks_taxon_filter.wdl" as taxon_filter
import "tasks_assembly.wdl" as assembly

workflow assemble_denovo_bulk {
  
  Array[File] reads_unmapped_bam_files

  scatter(reads_unmapped_bam in reads_unmapped_bam_files) {
    call taxon_filter.filter_to_taxon {
      input:
        reads_unmapped_bam = reads_unmapped_bam
    }
    
    call assembly.assemble {
      input:
        reads_unmapped_bam = filter_to_taxon.taxfilt_bam
    }
  }
}