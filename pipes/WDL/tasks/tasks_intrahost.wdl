version 1.0

task isnvs_per_sample {
  input {
    File    mapped_bam
    File    assembly_fasta

    Int?    threads
    Int?    minReadsPerStrand
    Int?    maxBias

    Int?    machine_mem_gb
    String  docker="quay.io/broadinstitute/viral-phylo:2.1.4.0"

    String  sample_name = basename(basename(basename(mapped_bam, ".bam"), ".all"), ".mapped")
  }

  command {
    intrahost.py --version | tee VERSION
    echo ${sample_name} | tee SAMPLE_NAME
    intrahost.py vphaser_one_sample \
        ${mapped_bam} \
        ${assembly_fasta} \
        ${sample_name}.vphaser2.txt.gz \
        ${'--vphaserNumThreads' + threads} \
        --removeDoublyMappedReads \
        ${'--minReadsEach' + minReadsPerStrand} \
        ${'--maxBias' + maxBias}
  }

  output {
    File   isnvsFile        = "${sample_name}.vphaser2.txt.gz"
    String sample_name_out  = read_string("SAMPLE_NAME")
    String viralngs_version = read_string("VERSION")
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    dx_instance_type: "mem1_ssd1_v2_x8"
  }
}


task isnvs_vcf {
  input {
    Array[File]    vphaser2Calls
    Array[File]    perSegmentMultiAlignments
    File           reference_fasta

    Array[String]? snpEffRef
    Array[String]? sampleNames
    String?        emailAddress
    Boolean        naiveFilter=false
    Boolean        append_reference_to_input=false
    Boolean        outputAllPositions=true

    Int?           machine_mem_gb
    String         docker="quay.io/broadinstitute/viral-phylo:2.1.4.0"
  }

  Int num_input_fastas = length(perSegmentMultiAlignments)

  parameter_meta {
    vphaser2Calls:             { description: "vphaser output; ex. vphaser2.<sample>.txt.gz" }
    perSegmentMultiAlignments: { description: "aligned_##.fasta, where ## is segment number" }
    snpEffRef:                 { description: "list of accessions to build/find snpEff database" }
    sampleNames:               { description: "list of sample names" }
    append_reference_to_input: { description: "Concatenate contents of reference sequence file to input sequence file, useful in the case where the input is already in the reference coordinate space and the reference sequence is not in the same file. The reference must contain exactly one segment since we cannot necessarily know without alignment which of the per-segment inputs each of several reference corresponds to." }
    emailAddress:              { description: "email address passed to NCBI if we need to download reference sequences" }
  }

  command {
    set -ex -o pipefail

    intrahost.py --version | tee VERSION

    if [[ "${append_reference_to_input}" == "true" ]]; then
      if [[ "${num_input_fastas}" == "1" ]]; then
        if [[ "$(grep '>' ${reference_fasta} | wc -l | tr -d ' ')" == "1" ]]; then
          cat ${sep=' ' perSegmentMultiAlignments} reference_fasta > ref_and_data_aligned_to_ref.fasta
          input_alignments="ref_and_data_aligned_to_ref.fasta"
        else
          echo "Reference fasta has >1 sequence, target for concatenation is unclear."
          exit 1
        fi
      else
        echo ">1 input fasta, unclear to which reference fasta data should be concatenated"
        exit 1
      fi
    else
      input_alignments="${sep=' ' perSegmentMultiAlignments}"
    fi

    SAMPLES="${sep=' ' sampleNames}"
    if [ -n "$SAMPLES" ]; then SAMPLES="--samples $SAMPLES"; fi

    providedSnpRefAccessions="${sep=' ' snpEffRef}"
    if [ -n "$providedSnpRefAccessions" ]; then 
      snpRefAccessions="$providedSnpRefAccessions";
    else
      snpRefAccessions="$(python -c "from Bio import SeqIO; print(' '.join(list(s.id for s in SeqIO.parse('${reference_fasta}', 'fasta'))))")"
    fi

    echo "snpRefAccessions: $snpRefAccessions"

    intrahost.py merge_to_vcf \
        ${reference_fasta} \
        isnvs.vcf.gz \
        $SAMPLES \
        --isnvs ${sep=' ' vphaser2Calls} \
        --alignments $input_alignments \
        --strip_chr_version \
        ${true="--naive_filter" false="" naiveFilter} \
        ${true="--output_all_positions" false="" outputAllPositions} \
        --parse_accession
        
    interhost.py snpEff \
        isnvs.vcf.gz \
        $snpRefAccessions \
        isnvs.annot.vcf.gz \
        ${'--emailAddress=' + emailAddress}

    intrahost.py iSNV_table \
        isnvs.annot.vcf.gz \
        isnvs.annot.txt.gz
  }

  output {
    File        isnvs_vcf           = "isnvs.vcf.gz"
    File        isnvs_vcf_idx       = "isnvs.vcf.gz.tbi"
    File        isnvs_annot_vcf     = "isnvs.annot.vcf.gz"
    File        isnvs_annot_vcf_idx = "isnvs.annot.vcf.gz.tbi"
    File        isnvs_annot_txt     = "isnvs.annot.txt.gz"
    String      viralngs_version    = read_string("VERSION")
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 4]) + " GB"
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}

task annotate_vcf_snpeff {
  input {
    File           in_vcf
    File           ref_fasta

    Array[String]? snpEffRef
    String?        emailAddress

    Int?           machine_mem_gb
    String         docker="quay.io/broadinstitute/viral-phylo:2.1.4.0"

    String         output_basename = basename(basename(in_vcf, ".gz"), ".vcf")
  }

  parameter_meta {
    in_vcf:             { description: "input VCF to annotate with snpEff", patterns: ["*.vcf","*.vcf.gz"]}
    ref_fasta:          { description: "The sequence containing the accession to use for annotation; only used if snpEffRef is not provided.", patterns: ["*.fasta","*.fa"] }
    snpEffRef:          { description: "list of accessions to build/find snpEff database. If this is not provided, the ID from the reference fasta will be used (it must be a GenBank accession)" }
    emailAddress:       { description: "email address passed to NCBI if we need to download reference sequences" }
  }

  command {
    set -ex -o pipefail

    intrahost.py --version | tee VERSION

    providedSnpRefAccessions="${sep=' ' snpEffRef}"
    if [ -n "$providedSnpRefAccessions" ]; then 
      snpRefAccessions="$providedSnpRefAccessions";
    else
      snpRefAccessions="$(python -c "from Bio import SeqIO; print(' '.join(list(s.id for s in SeqIO.parse('${ref_fasta}', 'fasta'))))")"
    fi
    echo "snpRefAccessions: $snpRefAccessions"

    if (file "${in_vcf}" | grep -q "gzip" ) ; then
      echo "${in_vcf} is already compressed"
    else
      echo "${in_vcf} is not compressed; gzipping..."
      bgzip "${in_vcf}"
    fi
    echo "Creating vcf index"
    tabix -p vcf "${in_vcf}"
        
    interhost.py snpEff \
        "${in_vcf}" \
        $snpRefAccessions \
        "${output_basename}.annot.vcf.gz" \
        ${'--emailAddress=' + emailAddress}

    intrahost.py iSNV_table \
        "${output_basename}.annot.vcf.gz" \
        "${output_basename}.annot.txt.gz"

    tabix -p vcf "${output_basename}.annot.vcf.gz"
  }

  output {
    File        annot_vcf_gz      = "${output_basename}.annot.vcf.gz"
    File        annot_vcf_gz_tbi  = "${output_basename}.annot.vcf.gz.tbi"
    File        annot_txt_gz      = "${output_basename}.annot.txt.gz"
    String      viralngs_version  = read_string("VERSION")
  }
  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 4]) + " GB"
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}
