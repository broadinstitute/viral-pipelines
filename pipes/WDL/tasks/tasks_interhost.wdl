version 1.0

task multi_align_mafft_ref {
  input {
    File           reference_fasta
    Array[File]+   assemblies_fasta # fasta files, one per sample, multiple chrs per file okay
    String         fasta_basename = basename(reference_fasta, '.fasta')
    Int?           mafft_maxIters
    Float?         mafft_ep
    Float?         mafft_gapOpeningPenalty

    Int?           machine_mem_gb
    String         docker="quay.io/broadinstitute/viral-phylo"
  }

  command {
    interhost.py --version | tee VERSION
    interhost.py multichr_mafft \
      ${reference_fasta} ${sep=' ' assemblies_fasta} \
      . \
      ${'--ep=' + mafft_ep} \
      ${'--gapOpeningPenalty=' + mafft_gapOpeningPenalty} \
      ${'--maxiters=' + mafft_maxIters} \
      --outFilePrefix align_mafft-${fasta_basename} \
      --preservecase \
      --localpair \
      --sampleNameListFile align_mafft-${fasta_basename}-sample_names.txt \
      --loglevel DEBUG
  }

  output {
    #File         sampleNamesFile    = "align_mafft-${fasta_basename}-sample_names.txt"
    Array[File]+  alignments_by_chr  = glob("align_mafft-${fasta_basename}*.fasta")
    String        viralngs_version   = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 28]) + " GB"
    cpu: 8
    dx_instance_type: "mem2_ssd1_v2_x8"
  }
}

task multi_align_mafft {
  input {
    Array[File]+   assemblies_fasta # fasta files, one per sample, multiple chrs per file okay
    String         out_prefix = "aligned"
    Int?           mafft_maxIters
    Float?         mafft_ep
    Float?         mafft_gapOpeningPenalty

    Int?           machine_mem_gb
    String         docker="quay.io/broadinstitute/viral-phylo"
  }

  command {
    interhost.py --version | tee VERSION
    interhost.py multichr_mafft \
      ${sep=' ' assemblies_fasta} \
      . \
      ${'--ep=' + mafft_ep} \
      ${'--gapOpeningPenalty=' + mafft_gapOpeningPenalty} \
      ${'--maxiters=' + mafft_maxIters} \
      --outFilePrefix ${out_prefix} \
      --preservecase \
      --localpair \
      --sampleNameListFile ${out_prefix}-sample_names.txt \
      --loglevel DEBUG
  }

  output {
    File        sampleNamesFile   = "${out_prefix}-sample_names.txt"
    Array[File] alignments_by_chr = glob("${out_prefix}*.fasta")
    String      viralngs_version  = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu: 8
    dx_instance_type: "mem1_ssd1_v2_x8"
  }
}

task index_ref {
  input {
    File     referenceGenome
    File?    novocraft_license

    Int?     machine_mem_gb
    String   docker="quay.io/broadinstitute/viral-core"
  }

  command {
    read_utils.py --version | tee VERSION
    read_utils.py novoindex \
    "${referenceGenome}" \
    ${"--NOVOALIGN_LICENSE_PATH=" + novocraft_license}
    
    read_utils.py index_fasta_samtools "${referenceGenome}"
    read_utils.py index_fasta_picard "${referenceGenome}"
  }

  output {
    File   referenceNix     = "*.nix"
    File   referenceFai     = "*.fasta.fai"
    File   referenceDict    = "*.dict"
    String viralngs_version = read_string("VERSION")
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 4]) + " GB"
  }
}

task trimal_clean_msa {
  input {
    File     in_aligned_fasta

    Int?     machine_mem_gb
    String   docker="quay.io/biocontainers/trimal:1.4.1--h6bb024c_3"

    String   input_basename = basename(basename(in_aligned_fasta, ".fasta"), ".fa")
  }

  command {
    trimal -fasta -automated1 -in "${in_aligned_fasta}" -out "${input_basename}_trimal_cleaned.fasta"
  }

  output {
    File   trimal_cleaned_fasta = "${input_basename}_trimal_cleaned.fasta"
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    dx_instance_type: "mem1_ssd1_v2_x8"
  }
}


