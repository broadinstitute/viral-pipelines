import "assemble_denovo.wdl" as denovo_assembly

workflow assemble_denovo_bulk {
  
  Array[File]+ reads_unmapped_bam_files
  File lastal_db_fasta
  Array[File]+ reference_genome_fasta
  trim_clip_db

  scatter(reads_unmapped_bam in reads_unmapped_bam_files) {
    call denovo_assembly.assemble_denovo {
      input:
        reads_unmapped_bam = reads_unmapped_bam,
        lastal_db_fasta = lastal_db_fasta,
        reference_genome_fasta = reference_genome_fasta,
        trim_clip_db = trim_clip_db
    }
  }
}