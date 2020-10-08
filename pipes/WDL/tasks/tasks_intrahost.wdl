version 1.0

task isnvs_per_sample {
  input {
    File    mapped_bam
    File    assembly_fasta

    Int?    threads
    Int?    minReadsPerStrand
    Int?    maxBias
    Boolean removeDoublyMappedReads=true

    Int?    machine_mem_gb
    String  docker="quay.io/broadinstitute/viral-phylo:2.1.4.0"

    String  sample_name = basename(basename(basename(mapped_bam, ".bam"), ".all"), ".mapped")
  }

  command {
    intrahost.py --version | tee VERSION
    intrahost.py vphaser_one_sample \
        ${mapped_bam} \
        ${assembly_fasta} \
        vphaser2.${sample_name}.txt.gz \
        ${'--vphaserNumThreads=' + threads} \
        ${true="--removeDoublyMappedReads" false="" removeDoublyMappedReads} \
        ${'--minReadsEach=' + minReadsPerStrand} \
        ${'--maxBias=' + maxBias}
  }

  output {
    File   isnvsFile        = "vphaser2.${sample_name}.txt.gz"
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

    Int?           machine_mem_gb
    String         docker="quay.io/broadinstitute/viral-phylo:2.1.4.0"
  }

  parameter_meta {
    vphaser2Calls:             { description: "vphaser output; ex. vphaser2.<sample>.txt.gz" }
    perSegmentMultiAlignments: { description: "aligned_##.fasta, where ## is segment number" }
    snpEffRef:                 { description: "list of accessions to build/find snpEff database" }
    sampleNames:               { description: "list of sample names" }
    emailAddress:              { description: "email address passed to NCBI if we need to download reference sequences" }
  }

  command {
    set -ex -o pipefail

    intrahost.py --version | tee VERSION

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
        --alignments ${sep=' ' perSegmentMultiAlignments} \
        --strip_chr_version \
        ${true="--naive_filter" false="" naiveFilter} \
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
