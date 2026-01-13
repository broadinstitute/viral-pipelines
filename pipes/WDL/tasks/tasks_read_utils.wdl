version 1.0

task max {
  input {
    Array[Int] list
    Int        default_empty = 0
  }
  Int disk_size = 10
  command <<<
    python3 << CODE
    inlist = '~{sep="*" list}'.split('*')
    print(str(max(map(int, [x for x in inlist if x]), default = ~{default_empty})))
    CODE
  >>>
  output {
    Int out = read_int(stdout())
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task group_bams_by_sample {
  input {
    Array[File] bam_filepaths
  }
  Int disk_size = 100
  parameter_meta {
    bam_filepaths: {
      description: "all bam files",
      localization_optional: true,
      stream: true,
      patterns: ["*.bam"]
    }
  }
  command <<<
    python3 << CODE
    import os.path

    # WDL arrays to python arrays
    bam_uris = '~{sep="*" bam_filepaths}'.split('*')

    # lookup table files to dicts
    sample_to_bams = {}
    for bam in bam_uris:
      # filename must be <samplename>.l<xxx>.bam
      assert bam.endswith('.bam'), "filename does not end in .bam: {}".format(bam)
      bam_base = os.path.basename(bam)
      i = bam_base.index('.l')
      assert i>0, "filename does not contain a .l -- {}".format(bam)
      sample = bam_base[:i]
      sample_to_bams.setdefault(sample, [])
      sample_to_bams[sample].append(bam)

    # write outputs
    with open('grouped_bams', 'wt') as out_groups:
      with open('sample_names', 'wt') as out_samples:
        for sample in sorted(sample_to_bams.keys()):
          out_samples.write(sample+'\n')
          out_groups.write('\t'.join(sample_to_bams[sample])+'\n')
    CODE
  >>>
  output {
    Array[Array[File]+] grouped_bam_filepaths = read_tsv('grouped_bams')
    Array[String]       sample_names          = read_lines('sample_names')
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task get_bam_samplename {
  input {
    File    bam
    String  docker = "quay.io/broadinstitute/viral-core:2.5.21"
  }
  Int   disk_size = round(size(bam, "GB")) + 50
  command <<<
    set -e -o pipefail
    samtools view -H "~{bam}" | \
      perl -lane 'if (/^\@RG\t.*SM:(\S+)/) { print "$1" }' | \
      sort | uniq > SAMPLE_NAME
  >>>
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
  output {
    String sample_name = read_string("SAMPLE_NAME")
  }
}

task get_sample_meta {
  input {
    Array[File] samplesheets_extended

    String      docker = "quay.io/broadinstitute/viral-core:2.5.21"
  }
  Int disk_size = 50
  command <<<
    python3 << CODE
    import os.path
    import csv
    import json
    import util.file

    # WDL arrays to python arrays
    library_metadata = '~{sep="*" samplesheets_extended}'.split('*')

    # lookup table files to dicts
    meta = {}
    meta_cols = ('sample','amplicon_set','control')
    for col in meta_cols:
      meta[col] = {}
    for libfile in library_metadata:
      with open(libfile, 'rt') as inf:
        for row in csv.DictReader(inf, delimiter='\t'):
          sanitized = util.file.string_to_file_name(row['sample'])
          for col in meta_cols:
            meta[col].setdefault(sanitized, '')
            if row.get(col):
              meta[col][sanitized] = row[col]

    # write outputs
    for col in meta_cols:
      with open(col, 'wt') as outf:
        json.dump(meta[col], outf, indent=2)
    CODE
  >>>
  output {
    Map[String,String] original_names = read_json('sample')
    Map[String,String] amplicon_set   = read_json('amplicon_set')
    Map[String,String] control        = read_json('control')
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}


task merge_and_reheader_bams {
    meta {
      description: "Merge and/or reheader bam files using a mapping table. This task can modify read group tags in a BAM header field for single BAM files or as part of a BAM merge operation. The output is a single BAM file (given one or more input BAMs) and a three-column tab delimited text table that defines: the field, the old value, and the new value (e.g. LB, old_lib_name, new_lib_name or SM, old_sample_name, new_sample_name)"
    }

    input {
      Array[File]+ in_bams
      String?      sample_name
      File?        reheader_table
      String       out_basename = basename(in_bams[0], ".bam")

      String       docker = "quay.io/broadinstitute/viral-core:2.5.21"
      Int          disk_size = 750
      Int          machine_mem_gb = 8
    }
    
    command <<<
        set -ex -o pipefail

        read_utils.py --version | tee VERSION
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

        if [ ~{length(in_bams)} -gt 1 ]; then
            read_utils.py merge_bams ~{sep=' ' in_bams} merged.bam --JVMmemory="$mem_in_mb"m --loglevel DEBUG \
                --picardOptions COMPRESSION_LEVEL=3 MAX_RECORDS_IN_RAM=200000 VALIDATION_STRINGENCY=SILENT
        else
            echo "Skipping merge, only one input file"
            cp ~{sep=' ' in_bams} merged.bam
        fi    

        # remap all SM values to user specified value
        if [ -n "~{sample_name}" ]; then
          # create sample name remapping table based on existing sample names
          samtools view -H merged.bam | perl -n -e'/SM:(\S+)/ && print "SM\t$1\t'"~{sample_name}"'\n"' | sort | uniq >> reheader_table.txt
        fi

        # remap arbitrary headers using user specified table
        if [[ -f "~{reheader_table}" ]]; then
          cat "~{reheader_table}" >> reheader_table.txt
        fi

        # reheader bam file if requested
        if [ -s reheader_table.txt ]; then
          read_utils.py reheader_bam merged.bam reheader_table.txt "~{out_basename}.bam" --loglevel DEBUG
        else
          mv merged.bam "~{out_basename}.bam"
        fi

        # summary stats on merged output
        samtools view -c "~{out_basename}.bam" | tee read_count_merged
        samtools flagstat "~{out_basename}.bam" | tee "~{out_basename}.bam.flagstat.txt"
    >>>

    output {
        File   out_bam          = "~{out_basename}.bam"
        Int    read_count       = read_int("read_count_merged")
        File   flagstat         = "~{out_basename}.bam.flagstat.txt"
        String viralngs_version = read_string("VERSION")
    }

    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        cpu: 4
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd2_v2_x4"
        preemptible: 0
        maxRetries: 2
    }
}

task rmdup_ubam {
  meta {
    description: "Perform read deduplication on unaligned reads."
  }

  input {
    File    reads_unmapped_bam
    String  method = "mvicuna"

    Int     max_reads = 100000000
    Int?    machine_mem_gb
    String  docker = "quay.io/broadinstitute/viral-core:2.5.21"
  }

  # Memory autoscaling: M-Vicuna loads reads into memory for deduplication.
  # For files <= 2GB, use 8GB RAM (covers typical samples).
  # For larger files, scale at 2x file size: 8 + 2*(file_size - 2) GB, capped at 128GB.
  # The 2x multiplier accounts for mvicuna's in-memory data structures which can
  # exceed the compressed BAM size, especially for deeply sequenced samples.
  Float input_size_gb = size(reads_unmapped_bam, "GB")
  Int mem_auto_scaled = if input_size_gb <= 2.0 then 8 else (if (8 + 2*ceil(input_size_gb - 2.0)) > 128 then 128 else (8 + 2*ceil(input_size_gb - 2.0)))
  Int mem_gb = select_first([machine_mem_gb, mem_auto_scaled])

  Int disk_size = ceil((3 * input_size_gb + 50) / 375.0) * 375

  parameter_meta {
    reads_unmapped_bam: { description: "unaligned reads in BAM format", patterns: ["*.bam"] }
    method:             { description: "mvicuna or cdhit" }
    max_reads:          { description: "If the input has more than this many reads, downsample before deduplication to avoid memory issues. Set to 0 to disable." }
  }

  String reads_basename = basename(reads_unmapped_bam, ".bam")

  command <<<
    set -ex -o pipefail
    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)
    read_utils.py --version | tee VERSION

    # Count input reads
    INPUT_BAM="~{reads_unmapped_bam}"
    READ_COUNT=$(samtools view -c "$INPUT_BAM")
    echo "$READ_COUNT" > dedup_read_count_pre
    echo "Input read count: $READ_COUNT"

    # Downsample if exceeds max_reads threshold
    if [ ~{max_reads} -gt 0 ] && [ "$READ_COUNT" -gt ~{max_reads} ]; then
      echo "Downsampling from $READ_COUNT to ~{max_reads} reads before deduplication"
      read_utils.py downsample_bams \
        "$INPUT_BAM" \
        --outPath ./downsample_out \
        --readCount=~{max_reads} \
        --JVMmemory "$mem_in_mb"m
      INPUT_BAM=$(ls downsample_out/*.bam)
      echo "Downsampled BAM: $INPUT_BAM"
    fi

    read_utils.py rmdup_"~{method}"_bam \
      "$INPUT_BAM" \
      "~{reads_basename}".dedup.bam \
      --loglevel=DEBUG

    samtools view -c "~{reads_basename}.dedup.bam" | tee dedup_read_count_post
  >>>

  output {
    File   dedup_bam             = "~{reads_basename}.dedup.bam"
    Int    dedup_read_count_pre  = read_int("dedup_read_count_pre")
    Int    dedup_read_count_post = read_int("dedup_read_count_post")
    String viralngs_version      = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu:    2
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem2_ssd1_v2_x2"
    maxRetries: 2
  }
}

task bbnorm_bam {
  meta {
    description: "Normalize read depth in a BAM file using BBNorm (kmer-based normalization). This uses BBMap's bbnorm tool to normalize coverage depth by downsampling over-represented kmers."
  }

  input {
    File    reads_bam

    Int     target = 1000
    Int?    kmer_length
    Int?    passes
    Int?    min_input_reads
    Int?    max_output_reads

    Int?    machine_mem_gb
    Int?    cpu
    String  docker = "quay.io/broadinstitute/viral-core:2.5.21"
  }

  # Memory autoscaling: BBNorm uses Java and loads kmer data structures into memory.
  # Scale from 8GB to 32GB based on input size. The 8-32GB range is cost-optimal on GCP
  # (1-4 GB/core with 8 CPUs). Observed peak usage was ~24GB for 36GB input files.
  Float input_size_gb = size(reads_bam, "GB")
  Int mem_auto_scaled = if input_size_gb <= 2.0 then 8 else (if (8 + ceil(input_size_gb - 2.0)) > 32 then 32 else (8 + ceil(input_size_gb - 2.0)))
  Int mem_gb = select_first([machine_mem_gb, mem_auto_scaled])

  # Passes autoscaling: For large inputs (>= 15GB), use 1 pass to optimize runtime over accuracy
  # For smaller inputs, use 2 passes for more accurate normalization.
  # Caller can override by specifying passes explicitly.
  Int passes_auto = if input_size_gb >= 15.0 then 1 else 2
  Int passes_to_use = select_first([passes, passes_auto])

  # Minimum 8 CPUs to maximize Google Cloud network bandwidth for file localization
  Int cpu_count = select_first([cpu, 8])

  # Disk: BAM->FASTQ expansion (~5-10x), BBNorm output + temp files (~10-15x more)
  # Note: GCP local SSDs must be allocated in pairs (2, 4, 8, 16, 24 × 375GB), so we round to 750GB multiples.
  Int disk_size = ceil((25 * input_size_gb + 50) / 750.0) * 750

  parameter_meta {
    reads_bam: {
      description: "Input reads in BAM format (may be aligned or unaligned, paired or unpaired).",
      patterns: ["*.bam"]
    }
    target: {
      description: "BBNorm target normalization depth. Reads are downsampled to achieve approximately this coverage depth. (default: 10000)"
    }
    kmer_length: {
      description: "Kmer length for BBNorm analysis. Longer kmers are more specific but require more memory. (default: bbnorm default of 31)"
    }
    passes: {
      description: "Number of normalization passes. More passes give more accurate normalization but take longer. (default: 2 for inputs < 15GB, 1 for inputs >= 15GB to optimize runtime)"
    }
    min_input_reads: {
      description: "If input has fewer than this many reads, skip normalization and copy input to output unchanged."
    }
    max_output_reads: {
      description: "If the normalized output would have more than this many read pairs, randomly downsample to this count. For paired-end data, the actual output read count will be 2x this value."
    }
  }

  String reads_basename = basename(reads_bam, ".bam")

  command <<<
    set -ex -o pipefail
    read_utils.py --version | tee VERSION

    # Count input reads first
    samtools view -c "~{reads_bam}" | tee bbnorm_read_count_pre

    # Calculate memory for BBNorm (85% of available, more accurate than bbnorm's auto-detect)
    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 85)

    read_utils.py rmdup_bbnorm_bam \
      "~{reads_bam}" \
      "~{reads_basename}.bbnorm.bam" \
      --target=~{target} \
      --passes=~{passes_to_use} \
      --memory="$mem_in_mb"m \
      ~{'--kmerLength=' + kmer_length} \
      ~{'--minInputReads=' + min_input_reads} \
      ~{'--maxOutputReads=' + max_output_reads} \
      --loglevel=DEBUG

    samtools view -c "~{reads_basename}.bbnorm.bam" | tee bbnorm_read_count_post
  >>>

  output {
    File   bbnorm_bam             = "~{reads_basename}.bbnorm.bam"
    Int    bbnorm_read_count_pre  = read_int("bbnorm_read_count_pre")
    Int    bbnorm_read_count_post = read_int("bbnorm_read_count_post")
    String viralngs_version       = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu:    cpu_count
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem2_ssd1_v2_x8"
    maxRetries: 2
  }
}

task downsample_bams {
  meta {
    description: "Downsample reads in a BAM file randomly subsampling to a target read count. Read deduplication can occur either before or after random subsampling, or not at all (default: not at all)."
    volatile: true
  }

  input {
    Array[File]+ reads_bam
    Int?         readCount
    Boolean      deduplicateBefore = false
    Boolean      deduplicateAfter = false

    Int?         machine_mem_gb
    String       docker = "quay.io/broadinstitute/viral-core:2.5.21"
  }

  Int disk_size = 750

  command <<<
    set -ex -o pipefail

    # find 90% memory
    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

    if [[ "~{deduplicateBefore}" == "true" ]]; then
      DEDUP_OPTION="--deduplicateBefore"
    elif [[ "~{deduplicateAfter}" == "true" ]]; then
      DEDUP_OPTION="--deduplicateAfter"
    fi

    if [[ "~{deduplicateBefore}" == "true" && "~{deduplicateAfter}" == "true" ]]; then
      echo "deduplicateBefore and deduplicateAfter are mutually exclusive. Only one can be used."
      exit 1
    fi

    read_utils.py --version | tee VERSION

    read_utils.py downsample_bams \
        ~{sep=' ' reads_bam} \
        --outPath ./output \
        ~{'--readCount=' + readCount} \
        $DEDUP_OPTION \
        --JVMmemory "$mem_in_mb"m
  >>>

  output {
    Array[File] downsampled_bam  = glob("output/*.downsampled-*.bam")
    String      viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 8]) + " GB"
    cpu:    4
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x4"
    maxRetries: 2
  }
}

task FastqToUBAM {
  meta {
    description: "Converts FASTQ (paired or single) to uBAM and adds read group information."
  }
  input {
    File    fastq_1
    File?   fastq_2
    String  sample_name
    String  library_name
    String? readgroup_name
    String? platform_unit
    String? run_date
    String  platform_name
    String? sequencing_center
    String? additional_picard_options

    Int     cpus = 2
    Int     mem_gb = 4
    Int     disk_size = 750
    String  docker = "quay.io/broadinstitute/viral-core:2.5.21"
  }
  parameter_meta {
    fastq_1: { description: "Unaligned read1 file in fastq format", patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"] }
    fastq_2: { description: "Unaligned read2 file in fastq format. This should be empty for single-end read conversion and required for paired-end reads. If provided, it must match fastq_1 in length and order.", patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"] }
    sample_name: { description: "Sample name. This is required and will populate the 'SM' read group value and will be used as the output filename (must be filename-friendly)." }
    library_name: { description: "Library name. This is required and will populate the 'LB' read group value. SM & LB combinations must be identical for any sequencing reads generated from the same sequencing library, and must be distinct for any reads generated from different libraries." }
    platform_name: { description: "Sequencing platform. This is required and will populate the 'PL' read group value. Must be one of CAPILLARY, DNBSEQ, HELICOS, ILLUMINA, IONTORRENT, LS454, ONT, PACBIO, or SOLID." }
    additional_picard_options: { description: "A string containing additional options to pass to picard FastqToSam beyond those made explicitly available as inputs to this task. For valid values, see: https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam" }
  }
  command <<<
      set -ex -o pipefail

      # find 90% memory
      mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

      read_utils.py --version | tee VERSION

      if [[ ! "~{platform_name}" =~ ^(CAPILLARY|DNBSEQ|ELEMENT|HELICOS|ILLUMINA|IONTORRENT|LS454|ONT|PACBIO|SINGULAR|SOLID|ULTIMA)$ ]]; then
        exit 1
      fi

      picard -Xmx"$mem_in_mb"m \
        FastqToSam \
        FASTQ="~{fastq_1}" \
        ~{"FASTQ2=" + fastq_2} \
        SAMPLE_NAME="~{sample_name}" \
        LIBRARY_NAME="~{library_name}" \
        OUTPUT="~{sample_name}".bam \
        ~{"READ_GROUP_NAME=" + readgroup_name} \
        ~{"PLATFORM_UNIT=" + platform_unit} \
        ~{"RUN_DATE=" + run_date} \
        ~{"PLATFORM=" + platform_name} \
        ~{"SEQUENCING_CENTER=" + sequencing_center} ~{additional_picard_options}
  >>>
  runtime {
    docker: docker
    cpu: cpus
    memory: mem_gb + " GB"
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
  output {
    File unmapped_bam = "~{sample_name}.bam"
  }
}

task read_depths {
  input {
    File      aligned_bam

    String    out_basename = basename(aligned_bam, '.bam')
    String    docker = "quay.io/broadinstitute/viral-core:2.5.21"
  }
  Int disk_size = 200
  command <<<
    set -e -o pipefail

    samtools depth "~{aligned_bam}" > "~{out_basename}.read_depths.txt"
  >>>

  output {
    File read_depths = "~{out_basename}.read_depths.txt"
  }
  runtime {
    docker: docker
    cpu:    2
    memory: "3 GB"
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}
