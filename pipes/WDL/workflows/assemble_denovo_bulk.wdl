import "assemble_denovo.wdl" as denovo_assembly

workflow assemble_denovo_bulk {
  
  Array[File]+ reads_unmapped_bam_files

  scatter(reads_unmapped_bam in reads_unmapped_bam_files) {
    call denovo_assembly.assemble_denovo {
      input: reads_unmapped_bam = reads_unmapped_bam
    }
  }
}