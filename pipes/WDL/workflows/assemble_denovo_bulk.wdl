import "assemble_denovo.wdl" as assemble_denovo

workflow assemble_denovo_bulk {
  
  Array[File]+ reads_unmapped_bam_files
  File lastal_db_fasta
  File trim_clip_db
  Array[File]+ reference_genome_fasta
  File? novocraft_license

  scatter(reads_unmapped_bam in reads_unmapped_bam_files) {
    call assemble_denovo_bulk.assemble_denovo {
      input:
        reads_unmapped_bam = reads_unmapped_bam,
        lastal_db_fasta = lastal_db_fasta,
        trim_clip_db = trim_clip_db,
        reference_genome_fasta = reference_genome_fasta,
        novocraft_license = novocraft_license
    }
  }
}