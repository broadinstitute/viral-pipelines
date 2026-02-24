version 1.0

task polyphonia_detect_cross_contamination {
  input {
    Array[File]  lofreq_vcfs
    Array[File]  genome_fastas
    File         reference_fasta

    Int          min_readcount       = 0
    Float        min_maf             = 0
    Float        min_genome_coverage = 0.95
    Int          min_read_depth      = 100
    Array[File]? read_depths
    String?      masked_positions
    File?        masked_positions_file
    Int          max_mismatches      = 0
    Int          min_matches         = 3
    Float        min_matches_proportion = 1

    Array[File]? plate_maps
    Int?         plate_size                 = 96
    Int?         plate_columns
    Int?         plate_rows
    Boolean?     compare_direct_neighbors   = true
    Boolean?     compare_diagonal_neighbors = false
    Boolean?     compare_full_row           = false
    Boolean?     compare_full_column        = false
    Boolean?     compare_full_plate         = false
    
    Boolean?     print_all_isnvs     = true
    Boolean?     print_all           = false

    String       out_basename = "potential_cross_contamination"

    String       docker = "quay.io/broadinstitute/polyphonia:latest"
  }

  Int disk_size = 100

  parameter_meta {
    lofreq_vcfs:         { description: "VCF file(s) output by LoFreq or GATK; must use reference provided by reference_fasta" }
    genome_fastas:       { description: "Unaligned consensus genome or genomes" }
    reference_fasta:     { description: "Reference fasta file" }
    
    min_readcount:       { description: "Minimum minor allele readcount for position to be considered heterozygous" }
    min_maf:             { description: "Minimum minor allele frequency for position to be considered heterozygous" }
    min_read_depth:      { description: "Minimum read depth for a position to be used for comparison" }
    read_depths:         { description: "Read depth tables (tab-separated, no header: reference name, position, read depth); provide alongside vcf files or heterozygosity tables if min-depth>0" }
    min_genome_coverage: { description: "Minimum proportion genome covered for a sample to be included" }
    masked_positions:    { description: "1-indexed positions to mask (e.g., 1-10,50,55-70)" }
    masked_positions_file: { description: "1-indexed positions to mask, one per line" }
    min_matches:         { description: "Of positions at which the two consensus genomes differ, the minimum number of positions at which contamination is detected as a minor allele" }
    min_matches_proportion: { description: "Of positions at which the two consensus genomes differ, the minimum proportion of positions at which contamination is detected as a minor allele" }
    max_mismatches:      { description: "Maximum allowed bases in contaminating sample consensus not matching contaminated sample alleles" }
    
    plate_maps:          { description: "Optional plate map(s) (tab-separated, no header: sample name, plate position (e.g., A8))" }
    plate_size:          { description: "Standard plate size (6-well, 12-well, 24, 48, 96, 384, 1536, 3456)" }
    plate_columns:       { description: "Number columns in plate (e.g., 1, 2, 3, 4)" }
    plate_rows:          { description: "Number rows in plate (e.g., A, B, C, D)" }
    compare_direct_neighbors:  { description: "Compare direct plate neighbors (left, right, top, bottom)" }
    compare_diagonal_neighbors:{ description: "Compare diagonal plate neighbors (top-right, bottom-right, top-left, bottom-left)" }
    compare_full_row:    { description: "Compare samples in the same row (e.g., row A)" }
    compare_full_column: { description: "Compare samples in the same column (e.g., column 8)" }
    compare_full_plate:  { description: "Compare all samples in the same plate map" }
    
    print_all_isnvs:     { description: "Include all threshold-passing samples in iSNVs visualizations, including samples without plate neighbors" }
    print_all:           { description: "Output outcomes of all comparisons (all comparisons are marked as potential cross-contamination)" }
  }

  command <<<
    set -e -o pipefail

    # commented out until polyphonia can report its own version
    #polyphonia --version | tee POLYPHONIA_VERSION

    mkdir -p figs
    
    # prep plate maps input, if specified
    PLATE_MAPS_INPUT=""
    if [ -n "${sep=' ' select_first([plate_maps, []])}" ]; then
      PLATE_MAPS_INPUT="--plate-map ~{sep=' ' select_first([plate_maps, []])}"
    fi

    polyphonia cross_contamination \
      --ref ~{reference_fasta} \
      --vcf ~{sep=' ' lofreq_vcfs} \
      --consensus ~{sep=' ' genome_fastas} \
      --read-depths ~{sep=' ' select_first([read_depths, []])} \
      ~{'--min-covered ' + min_genome_coverage} \
      ~{'--min-depth ' + min_read_depth} \
      ~{'--min-readcount ' + min_readcount} \
      ~{'--max-mismatches ' + max_mismatches} \
      ~{'--min-matches ' + min_matches} \
      ~{'--min-matches-proportion ' + min_matches_proportion} \
      ~{'--min-maf ' + min_maf} \
      ~{'--masked-positions ' + masked_positions} \
      ~{'--masked-positions-file ' + masked_positions_file} \
      $PLATE_MAPS_INPUT \
      ~{'--plate-size ' + plate_size} \
      ~{'--plate-columns ' + plate_columns} \
      ~{'--plate-rows ' + plate_rows} \
      --compare-direct ~{if select_first([compare_direct_neighbors, true]) then "TRUE" else "FALSE"} \
      --compare-diagonal ~{if select_first([compare_diagonal_neighbors, false]) then "TRUE" else "FALSE"} \
      --compare-row ~{if select_first([compare_full_row, false]) then "TRUE" else "FALSE"} \
      --compare-column ~{if select_first([compare_full_column, false]) then "TRUE" else "FALSE"} \
      --compare-plate ~{if select_first([compare_full_plate, false]) then "TRUE" else "FALSE"} \
      --output ~{out_basename}.txt \
      --out-figures figs \
      --cores `nproc` \
      --verbose TRUE \
      --overwrite TRUE \
      --print-all-iSNVs ~{if select_first([print_all_isnvs, true]) then "TRUE" else "FALSE"} \
      --print-all ~{if select_first([print_all, false]) then "TRUE" else "FALSE"}
  >>>

  output {
    File        report             = "~{out_basename}.txt"
    Array[File] figures            = glob("figs/*")
    # commented out until polyphonia can report its own version
    #String      polyphonia_version = read_string("POLYPHONIA_VERSION")
  }
  runtime {
    docker: docker
    cpu:    4
    memory: "13 GB"
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}

task lofreq {
  input {
    File      aligned_bam
    File      reference_fasta

    String    out_basename = basename(aligned_bam, '.bam')
    String    docker = "ghcr.io/broadinstitute/viral-ngs:3.0.6-phylo"
  }
  Int disk_size = 200
  command <<<
    set -e -o pipefail

    lofreq version | grep version | sed 's/.* \(.*\)/\1/g' | tee LOFREQ_VERSION

    # make local copies because CWD is writeable but localization dir isn't always
    cp "~{aligned_bam}" aligned.bam
    python3<<CODE
    import shutil
    from viral_ngs.core import file as util_file
    with util_file.fastas_with_sanitized_ids("~{reference_fasta}", use_tmp=True) as sanitized_fastas:
        shutil.copyfile(sanitized_fastas[0], 'reference.fasta')
    CODE

    # samtools faidx fails if fasta is empty
    if [ $(grep -v '^>' reference.fasta | tr -d '\nNn' | wc -c) == "0" ]; then
      touch "~{out_basename}.vcf"
      exit 0
    fi

    # index for lofreq
    samtools faidx reference.fasta
    samtools index aligned.bam

    # lofreq
    lofreq call \
      -f reference.fasta \
      -o "~{out_basename}.vcf" \
      aligned.bam
  >>>

  output {
    File   report_vcf     = "~{out_basename}.vcf"
    String lofreq_version = read_string("LOFREQ_VERSION")
  }
  runtime {
    docker: docker
    cpu:    2
    memory: "3 GB"
    disks: "local-disk ~{disk_size} SSD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}

task isnvs_per_sample {
  input {
    File    mapped_bam
    File    assembly_fasta

    Int?    threads
    Int?    minReadsPerStrand
    Int?    maxBias
    Boolean removeDoublyMappedReads = true

    Int?    machine_mem_gb
    String  docker = "ghcr.io/broadinstitute/viral-ngs:3.0.6-phylo"

    String  sample_name = basename(basename(basename(mapped_bam, ".bam"), ".all"), ".mapped")
  }
  

  command <<<
    intrahost --version | tee VERSION
    intrahost vphaser_one_sample \
        ~{mapped_bam} \
        ~{assembly_fasta} \
        vphaser2.~{sample_name}.txt.gz \
        ~{'--vphaserNumThreads=' + threads} \
        ~{true="--removeDoublyMappedReads" false="" removeDoublyMappedReads} \
        ~{'--minReadsEach=' + minReadsPerStrand} \
        ~{'--maxBias=' + maxBias}
  >>>

  output {
    File   isnvsFile        = "vphaser2.~{sample_name}.txt.gz"
    String viralngs_version = read_string("VERSION")
  }
  runtime {
    docker: docker
    memory: "~{select_first([machine_mem_gb, 7])} GB"
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
    Boolean        naiveFilter = false

    Int?           machine_mem_gb
    String         docker = "ghcr.io/broadinstitute/viral-ngs:3.0.6-phylo"
  }

  parameter_meta {
    vphaser2Calls:             { description: "vphaser output; ex. vphaser2.<sample>.txt.gz" }
    perSegmentMultiAlignments: { description: "aligned_##.fasta, where ## is segment number" }
    snpEffRef:                 { description: "list of accessions to build/find snpEff database" }
    sampleNames:               { description: "list of sample names" }
    emailAddress:              { description: "email address passed to NCBI if we need to download reference sequences" }
  }

  command <<<
    set -ex -o pipefail

    intrahost --version | tee VERSION

    SAMPLES="~{sep=' ' select_first([sampleNames, []])}"
    if [ -n "$SAMPLES" ]; then SAMPLES="--samples $SAMPLES"; fi

    providedSnpRefAccessions="~{sep=' ' select_first([snpEffRef, []])}"
    if [ -n "$providedSnpRefAccessions" ]; then
      snpRefAccessions="$providedSnpRefAccessions";
    else
      snpRefAccessions="$(python -c "from Bio import SeqIO; print(' '.join(list(s.id for s in SeqIO.parse('~{reference_fasta}', 'fasta'))))")"
    fi

    echo "snpRefAccessions: $snpRefAccessions"

    intrahost merge_to_vcf \
        ~{reference_fasta} \
        isnvs.vcf.gz \
        $SAMPLES \
        --isnvs ~{sep=' ' vphaser2Calls} \
        --alignments ~{sep=' ' perSegmentMultiAlignments} \
        --strip_chr_version \
        ~{true="--naive_filter" false="" naiveFilter} \
        --parse_accession

    interhost snpEff \
        isnvs.vcf.gz \
        $snpRefAccessions \
        isnvs.annot.vcf.gz \
        ~{'--emailAddress=' + emailAddress}

    intrahost iSNV_table \
        isnvs.annot.vcf.gz \
        isnvs.annot.txt.gz
  >>>

  output {
    File   isnvs_vcf           = "isnvs.vcf.gz"
    File   isnvs_vcf_idx       = "isnvs.vcf.gz.tbi"
    File   isnvs_annot_vcf     = "isnvs.annot.vcf.gz"
    File   isnvs_annot_vcf_idx = "isnvs.annot.vcf.gz.tbi"
    File   isnvs_annot_txt     = "isnvs.annot.txt.gz"
    String viralngs_version    = read_string("VERSION")
  }
  runtime {
    docker: docker
    memory: "~{select_first([machine_mem_gb, 4])} GB"
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
    String         docker = "ghcr.io/broadinstitute/viral-ngs:3.0.6-phylo"

    String         output_basename = basename(basename(in_vcf, ".gz"), ".vcf")
  }

  Int disk_size = 375


  parameter_meta {
    in_vcf:             { description: "input VCF to annotate with snpEff", patterns: ["*.vcf","*.vcf.gz"]}
    ref_fasta:          { description: "The sequence containing the accession to use for annotation; only used if snpEffRef is not provided.", patterns: ["*.fasta","*.fa"] }
    snpEffRef:          { description: "list of accessions to build/find snpEff database. If this is not provided, the ID from the reference fasta will be used (it must be a GenBank accession)" }
    emailAddress:       { description: "email address passed to NCBI if we need to download reference sequences" }
  }

  command <<<
    set -ex -o pipefail

    intrahost --version | tee VERSION

    providedSnpRefAccessions="~{sep=' ' select_first([snpEffRef, []])}"
    if [ -n "$providedSnpRefAccessions" ]; then
      snpRefAccessions="$providedSnpRefAccessions";
    else
      snpRefAccessions="$(python -c "from Bio import SeqIO; print(' '.join(list(s.id for s in SeqIO.parse('~{ref_fasta}', 'fasta'))))")"
    fi
    echo "snpRefAccessions: $snpRefAccessions"

    vcf_to_use=""
    if (file "~{in_vcf}" | grep -q "gzip" ) ; then
      echo "~{in_vcf} is already compressed"
      vcf_to_use="~{in_vcf}"
    else
      echo "~{in_vcf} is not compressed; gzipping..."
      bgzip "~{in_vcf}"
      vcf_to_use="~{in_vcf}.gz"
    fi

    # renames the seq id using the first sequence in the alignment
    ref_name=$(head -n1 "~{ref_fasta}" | sed -E 's/^>([^[:space:]]+).*$/\1/g')
    ref_name_no_version=$(head -n1 "~{ref_fasta}" | sed -E 's/^>([^[:space:]\.]+).*$/\1/g')
    # copy the input or just created gzipped vcf file
    cp "$vcf_to_use" "temp.vcf.gz"
    # ensure uncompressed
    bgzip -d "temp.vcf.gz"
    # rename chr field (first col) in vcf
    cat "temp.vcf" | sed "s/^1/$ref_name_no_version/" > "temp2.vcf"

    # output the vcf, removing the reference sequence if present as a sample name
    bgzip "temp2.vcf"
    tabix -p vcf "temp2.vcf.gz"
    bcftools index "temp2.vcf.gz"
    bcftools view -s "^$ref_name" "temp2.vcf.gz" > "temp3.vcf"
    rm "temp2.vcf.gz"
    vcf_to_use="temp3.vcf"

    # compress vcf
    bgzip -c "$vcf_to_use" > "$vcf_to_use.gz"
    vcf_to_use="$vcf_to_use.gz"

    # index vcf
    echo "Creating vcf index"
    bcftools index "$vcf_to_use"
    tabix -p vcf "$vcf_to_use"

    interhost snpEff \
        "$vcf_to_use" \
        $snpRefAccessions \
        "~{output_basename}.annot.vcf.gz" \
        ~{'--emailAddress=' + emailAddress}
  >>>

  output {
    File   annot_vcf_gz     = "~{output_basename}.annot.vcf.gz"
    File   annot_vcf_gz_tbi = "~{output_basename}.annot.vcf.gz.tbi"
    String viralngs_version = read_string("VERSION")
  }
  runtime {
    docker: docker
    memory: "~{select_first([machine_mem_gb, 4])} GB"
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}
