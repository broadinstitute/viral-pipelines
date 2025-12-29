version 1.0

task deplete_taxa {
  meta { description: "Runs a full human read depletion pipeline and removes PCR duplicates. Input database files (minimapDbs, bmtaggerDbs, blastDbs, bwaDbs) may be any combination of: .fasta, .fasta.gz, or tarred up indexed fastas (using the software's indexing method) as .tar.gz, .tar.bz2, .tar.lz4, or .tar.zst." }

  input {
    File         raw_reads_unmapped_bam
    Array[File]? minimapDbs
    Array[File]? bmtaggerDbs
    Array[File]? blastDbs
    Array[File]? bwaDbs
    Int?         query_chunk_size
    Boolean      clear_tags = false
    String       tags_to_clear_space_separated = "XT X0 X1 XA AM SM BQ CT XN OC OP"

    Int?         cpu
    Int?         machine_mem_gb
    String       docker = "quay.io/broadinstitute/viral-classify:2.5.16.0"
  }

  # Autoscale CPU based on input size: 8 CPUs for ~1M reads (0.15 GB), 96 CPUs for ~100M reads (15 GB)
  # Linear scaling: 8 + (input_GB / 15) * 88, capped at 96, rounded to nearest multiple of 4
  Float        input_bam_size_gb = size(raw_reads_unmapped_bam, "GB")
  Float        cpu_unclamped = 8.0 + (input_bam_size_gb / 15.0) * 88.0
  Int          cpu_actual = select_first([cpu, floor(((if cpu_unclamped > 96.0 then 96.0 else cpu_unclamped) + 2.0) / 4.0) * 4])
  # Memory scales with CPU at 2x ratio (default), or use override
  Int          machine_mem_gb_actual = select_first([machine_mem_gb, cpu_actual * 2])
  # Disable preemptible for large inputs (>6GB) to avoid restart delays on long-running jobs
  Int          preemptible_tries = if input_bam_size_gb > 6.0 then 0 else 1

  parameter_meta {
    raw_reads_unmapped_bam: { description: "unaligned reads in BAM format", patterns: ["*.bam"] }
    minimapDbs: {
       description: "Optional list of databases to use for minimap2-based depletion. Sequences in fasta format will be indexed on the fly, pre-indexed databases may be provided as tarballs.",
       patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    bmtaggerDbs: {
       description: "Optional list of databases to use for bmtagger-based depletion. Sequences in fasta format will be indexed on the fly, pre-bmtagger-indexed databases may be provided as tarballs.",
       patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    blastDbs: {
      description: "Optional list of databases to use for blastn-based depletion. Sequences in fasta format will be indexed on the fly, pre-blast-indexed databases may be provided as tarballs.",
      patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
    bwaDbs: {
      description: "Optional list of databases to use for bwa mem-based depletion. Sequences in fasta format will be indexed on the fly, pre-bwa-indexed databases may be provided as tarballs.",
      patterns: ["*.fasta", "*.fasta.gz", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
    }
  }

  String       bam_basename = basename(raw_reads_unmapped_bam, ".bam")
  Float        minimap_db_size = if defined(minimapDbs) then size(select_first([minimapDbs, []]), "GB") else 0
  Float        bmtagger_db_size = if defined(bmtaggerDbs) then size(select_first([bmtaggerDbs, []]), "GB") else 0
  Float        blast_db_size = if defined(blastDbs) then size(select_first([blastDbs, []]), "GB") else 0
  Float        bwa_db_size = if defined(bwaDbs) then size(select_first([bwaDbs, []]), "GB") else 0
  Float        total_db_size = minimap_db_size + bmtagger_db_size + blast_db_size + bwa_db_size
  Int          disk_size = ceil((10 * size(raw_reads_unmapped_bam, "GB") + 2 * total_db_size + 100) / 375.0) * 375

  command <<<
    set -ex -o pipefail
    taxon_filter.py --version | tee VERSION

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi

    # find memory thresholds
    mem_in_mb_50=$(/opt/viral-ngs/source/docker/calc_mem.py mb 50)
    mem_in_mb_75=$(/opt/viral-ngs/source/docker/calc_mem.py mb 75)

    # depletion db args
    DBS_MINIMAP="~{sep=' ' minimapDbs}"
    DBS_BMTAGGER="~{sep=' ' bmtaggerDbs}"
    DBS_BLAST="~{sep=' ' blastDbs}"
    DBS_BWA="~{sep=' ' bwaDbs}"
    if [ -n "$DBS_MINIMAP" ]; then DBS_MINIMAP="--minimapDbs $DBS_MINIMAP"; fi
    if [ -n "$DBS_BMTAGGER" ]; then DBS_BMTAGGER="--bmtaggerDbs $DBS_BMTAGGER"; fi
    if [ -n "$DBS_BLAST" ]; then DBS_BLAST="--blastDbs $DBS_BLAST"; fi
    if [ -n "$DBS_BWA" ]; then DBS_BWA="--bwaDbs $DBS_BWA"; fi
    
    if [[ "~{clear_tags}" == "true" ]]; then
      TAGS_TO_CLEAR="--clearTags"
      if [[ -n "~{tags_to_clear_space_separated}" ]]; then
        TAGS_TO_CLEAR="$TAGS_TO_CLEAR ~{'--tagsToClear=' + tags_to_clear_space_separated}"
      fi
    fi

    # run depletion
    taxon_filter.py deplete \
      "~{raw_reads_unmapped_bam}" \
      tmpfile.raw.bam \
      tmpfile.minimap.bam \
      tmpfile.bwa.bam \
      tmpfile.bmtagger_depleted.bam \
      "~{bam_basename}.cleaned.bam" \
      $DBS_MINIMAP $DBS_BMTAGGER $DBS_BLAST $DBS_BWA \
      ~{'--chunkSize=' + query_chunk_size} \
      $TAGS_TO_CLEAR \
      --JVMmemory="$mem_in_mb_75"m \
      --srprismMemory=$mem_in_mb_75 \
      --loglevel=DEBUG
    rm -f tmpfile.raw.bam tmpfile.minimap.bam tmpfile.bwa.bam tmpfile.bmtagger_depleted.bam

    samtools view -c "~{raw_reads_unmapped_bam}" | tee depletion_read_count_pre
    samtools view -c "~{bam_basename}.cleaned.bam" | tee depletion_read_count_post

    cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
    cat /proc/loadavg | cut -f 3 -d ' ' > LOAD_15M
    set +o pipefail
    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi; } > MEM_BYTES
  >>>

  output {
    File   cleaned_bam               = "~{bam_basename}.cleaned.bam"
    Int    depletion_read_count_pre  = read_int("depletion_read_count_pre")
    Int    depletion_read_count_post = read_int("depletion_read_count_post")
    Int    max_ram_gb                = ceil(read_float("MEM_BYTES")/1000000000)
    Int    runtime_sec               = ceil(read_float("UPTIME_SEC"))
    Int    cpu_load_15min            = ceil(read_float("LOAD_15M"))
    String viralngs_version          = read_string("VERSION")
  }
  runtime {
    docker: docker
    memory: machine_mem_gb_actual + " GB"
    cpu: cpu_actual
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x8"
    preemptible: preemptible_tries
    maxRetries: 1
  }
}

task filter_to_taxon {
  meta { description: "This step reduces the read set to a specific taxon (usually the genus level or greater for the virus of interest)" }

  input {
    File     reads_unmapped_bam
    File     lastal_db_fasta
    Boolean  error_on_reads_in_neg_control = false
    Int      negative_control_reads_threshold = 0
    String   neg_control_prefixes_space_separated = "neg water NTC"

    Int      machine_mem_gb = 15
    String   docker = "quay.io/broadinstitute/viral-classify:2.5.16.0"
  }

  # do this in two steps in case the input doesn't actually have "cleaned" in the name
  String   bam_basename = basename(basename(reads_unmapped_bam, ".bam"), ".cleaned")
  Int disk_size = ceil((6 * size(reads_unmapped_bam, "GB") + 2 * size(lastal_db_fasta, "GB") + 100) / 375.0) * 375

  command <<<
    set -ex -o pipefail
    taxon_filter.py --version | tee VERSION

    # find 90% memory
    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

    if [[ "~{error_on_reads_in_neg_control}" == "true" ]]; then
      ERROR_ON_NEG_CONTROL_ARGS="--errorOnReadsInNegControl"
      if [[ -n "~{negative_control_reads_threshold}" ]]; then
        ERROR_ON_NEG_CONTROL_ARGS="$ERROR_ON_NEG_CONTROL_ARGS ~{'--negativeControlReadsThreshold=' + negative_control_reads_threshold}"
      fi
      if [[ -n "~{neg_control_prefixes_space_separated}" ]]; then
        ERROR_ON_NEG_CONTROL_ARGS="$ERROR_ON_NEG_CONTROL_ARGS ~{'--negControlPrefixes=' + neg_control_prefixes_space_separated}"
      fi      
    fi

    taxon_filter.py filter_lastal_bam \
      "~{reads_unmapped_bam}" \
      "~{lastal_db_fasta}" \
      "~{bam_basename}.taxfilt.bam" \
      $ERROR_ON_NEG_CONTROL_ARGS \
      --JVMmemory="$mem_in_mb"m \
      --loglevel=DEBUG

    samtools view -c "~{bam_basename}.taxfilt.bam" | tee filter_read_count_post
  >>>

  output {
    File   taxfilt_bam            = "~{bam_basename}.taxfilt.bam"
    Int    filter_read_count_post = read_int("filter_read_count_post")
    String viralngs_version       = read_string("VERSION")
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    cpu: 16
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x8"
    maxRetries: 2
  }
}

task build_lastal_db {
  input {
    File   sequences_fasta

    Int    machine_mem_gb = 7
    String docker = "quay.io/broadinstitute/viral-classify:2.5.16.0"
  }

  String db_name = basename(sequences_fasta, ".fasta")
  Int    disk_size = 375

  command <<<
    set -ex -o pipefail
    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    taxon_filter.py --version | tee VERSION
    taxon_filter.py lastal_build_db "~{sequences_fasta}" ./ --loglevel=DEBUG
    tar -c ~{db_name}* | lz4 -9 > "~{db_name}.tar.lz4"
  >>>

  output {
    File   lastal_db        = "~{db_name}.tar.lz4"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x4"
    maxRetries: 2
  }
}

task merge_one_per_sample {
  input {
    String       out_bam_basename
    Array[File]+ inputBams
    Boolean      rmdup = false

    Int          machine_mem_gb = 7
    String       docker = "quay.io/broadinstitute/viral-core:2.5.17"
  }

  Int disk_size = 750

  command <<<
    set -ex -o pipefail
    read_utils.py --version | tee VERSION

    # find 90% memory
    mem_in_mb=~(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

    read_utils.py merge_bams \
      "~{sep=' ' inputBams}" \
      "~{out_bam_basename}.bam" \
      --picardOptions SORT_ORDER=queryname \
      --JVMmemory "$mem_in_mb"m \
      --loglevel=DEBUG

    if [[ "~{rmdup}" == "true" ]]; then
      mv "~{out_bam_basename}.bam" tmp.bam
      read_utils.py rmdup_mvicuna_bam \
        tmp.bam \
        "~{out_bam_basename}.bam" \
        --JVMmemory "$mem_in_mb"m \
        --loglevel=DEBUG
    fi
  >>>

  output {
    File   mergedBam        = "~{out_bam_basename}.bam"
    String viralngs_version = read_string("VERSION")
  }

  runtime{
    memory: machine_mem_gb + " GB"
    cpu: 4
    docker: docker
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd2_v2_x4"
    maxRetries: 2
  }
}
