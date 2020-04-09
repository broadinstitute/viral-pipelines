
task merge_tarballs {
  Array[File]+  tar_chunks
  String        out_filename
  String?       docker="quay.io/broadinstitute/viral-core"

#  parameter_meta {
#    tar_chunks: "stream"
#  }

  command {
    set -ex -o pipefail

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi

    file_utils.py --version | tee VERSION

    file_utils.py merge_tarballs \
      ${out_filename} ${sep=' ' tar_chunks} \
      --loglevel=DEBUG
  }

  output {
    File    combined_tar      = "${out_filename}"
    String  viralngs_version  = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "7 GB"
    cpu: 16
    disks: "local-disk 1125 LOCAL"
    dx_instance_type: "mem1_ssd2_v2_x16"
    preemptible: 0
  }
}

task illumina_demux {

  File    flowcell_tgz
  Int?    lane=1
  File?   samplesheet
  File?   runinfo
  String? sequencingCenter

  String? flowcell
  Int?    minimumBaseQuality = 10
  Int?    maxMismatches = 0
  Int?    minMismatchDelta
  Int?    maxNoCalls
  String? readStructure
  Int?    minimumQuality
  Int?    threads
  String? runStartDate
  Int?    maxReadsInRamPerTile
  Int?    maxRecordsInRam
  Boolean? forceGC=true

  String? docker="quay.io/broadinstitute/viral-core"


#  parameter_meta {
#    flowcell_tgz: "stream"
#  }

  command {
    set -ex -o pipefail

    # find N% memory
    mem_in_mb=`/opt/viral-ngs/source/docker/calc_mem.py mb 85`

    if [ -d /mnt/tmp ]; then
      TMPDIR=/mnt/tmp
    fi
    FLOWCELL_DIR=$(mktemp -d)

    read_utils.py --version | tee VERSION

    read_utils.py extract_tarball \
      ${flowcell_tgz} $FLOWCELL_DIR \
      --loglevel=DEBUG

    # if we are overriding the RunInfo file, use the path of the file provided. Otherwise find the file
    if [ -n "${runinfo}" ]; then
      RUNINFO_FILE="${runinfo}"
    else
      # full RunInfo.xml path
      RUNINFO_FILE="$(find $FLOWCELL_DIR -type f -maxdepth 3 -name RunInfo.xml | head -n 1)"
    fi
    
    # Parse the lane count & run ID from RunInfo.xml file
    lane_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@LaneCount)" $RUNINFO_FILE)
    if [ -z "$lane_count" ]; then
        echo "Could not parse LaneCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
    fi

    surface_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@SurfaceCount)" $RUNINFO_FILE)
    if [ -z "$surface_count" ]; then
        echo "Could not parse SurfaceCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
    fi

    swath_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@SwathCount)" $RUNINFO_FILE)
    if [ -z "$swath_count" ]; then
        echo "Could not parse SwathCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
    fi

    tile_count=$(xmllint --xpath "string(//Run/FlowcellLayout/@TileCount)" $RUNINFO_FILE)
    if [ -z "$tile_count" ]; then
        echo "Could not parse TileCount from RunInfo.xml. Please check RunInfo.xml is properly formatted"
    fi

    # total data size more roughly tracks total tile count
    total_tile_count=$((lane_count*surface_count*swath_count*tile_count))

    demux_threads="$(nproc --all)"
    if [ "$total_tile_count" -le 50 ]; then
        echo "Detected $total_tile_count tiles, interpreting as MiSeq run."
    elif [ "$total_tile_count" -le 150 ]; then
        echo "Detected $total_tile_count tiles, interpreting as HiSeq2k run."
    elif [ "$total_tile_count" -le 288 ]; then
        # increase the number of reads in ram per-tile for NextSeq, since the tiles are larger
        # without this setting, reads will spill to disk and may read the limit
        # on the number of files that can be opened
        max_reads_in_ram_per_tile=1500000
        max_records_in_ram=2000000
        echo "Detected $total_tile_count tiles, interpreting as NextSeq (mid-output) run."
    elif [ "$total_tile_count" -le 864 ]; then
        # increase the number of reads in ram per-tile for NextSeq, since the tiles are larger
        # without this setting, reads will spill to disk and may read the limit
        # on the number of files that can be opened
        max_reads_in_ram_per_tile=1500000 # reduce the number of reads per tile since the NovaSeq has so many
        max_records_in_ram=2500000
        echo "Detected $total_tile_count tiles, interpreting as NextSeq (high-output) run."
    elif [ "$total_tile_count" -le 896 ]; then
        echo "Detected $total_tile_count tiles, interpreting as HiSeq4k run."
    elif [ "$total_tile_count" -le 1408 ]; then
        demux_threads=32 # with NovaSeq-size output, OOM errors can sporadically occur with higher thread counts
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 80)
        max_reads_in_ram_per_tile=1200000 # reduce the number of reads per tile since the NovaSeq has so many
        max_records_in_ram=3000000
        echo "Detected $total_tile_count tiles, interpreting as NovaSeq run."
        echo "  **Note: Q20 threshold used since NovaSeq with RTA3 writes only four Q-score values: 2, 12, 23, and 37.**"
        echo "    See: https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/novaseq-hiseq-q30-app-note-770-2017-010.pdf"
    elif [ "$total_tile_count" -gt 1408 ]; then
        demux_threads=32 # with NovaSeq-size output, OOM errors can sporadically occur with higher thread counts
        max_reads_in_ram_per_tile=1200000 # reduce the number of reads per tile since the NovaSeq has so many
        max_records_in_ram=3000000
        echo "Tile count: $total_tile_count tiles (unknown instrument type)."
    fi

    # use the passed-in (or default) WDL value first, then fall back to the auto-scaled value
    # if the result of this is null (nothing is passed in, no autoscaled value, no param is passed to the command)
    if [ -n "${minimumBaseQuality}" ]; then demux_min_base_quality="${minimumBaseQuality}"; else demux_min_base_quality="$demux_min_base_quality"; fi
    if [ -n "$demux_min_base_quality" ]; then demux_min_base_quality="--minimum_base_quality=$demux_min_base_quality";fi
    
    if [ -n "${threads}" ]; then demux_threads="${threads}"; else demux_threads="$demux_threads"; fi
    if [ -n "$demux_threads" ]; then demux_threads="--threads=$demux_threads"; fi
    

    if [ -n "${maxReadsInRamPerTile}" ]; then max_reads_in_ram_per_tile="${maxReadsInRamPerTile}"; else max_reads_in_ram_per_tile="$max_reads_in_ram_per_tile"; fi
    if [ -n "$max_reads_in_ram_per_tile" ]; then max_reads_in_ram_per_tile="--max_reads_in_ram_per_tile=$max_reads_in_ram_per_tile"; fi
    
    if [ -n "${maxRecordsInRam}" ]; then max_records_in_ram="${maxRecordsInRam}"; else max_records_in_ram="$max_records_in_ram"; fi
    if [ -n "$max_records_in_ram" ]; then max_records_in_ram="--max_records_in_ram=$max_records_in_ram"; fi

    # note that we are intentionally setting --threads to about 2x the core
    # count. seems to still provide speed benefit (over 1x) when doing so.
    illumina.py illumina_demux \
      $FLOWCELL_DIR \
      ${lane} \
      . \
      ${'--sampleSheet=' + samplesheet} \
      ${'--runInfo=' + runinfo} \
      ${'--sequencing_center=' + sequencingCenter} \
      --outMetrics=metrics.txt \
      --commonBarcodes=barcodes.txt \
      ${'--flowcell=' + flowcell} \
      $demux_min_base_quality \
      ${'--max_mismatches=' + maxMismatches} \
      ${'--min_mismatch_delta=' + minMismatchDelta} \
      ${'--max_no_calls=' + maxNoCalls} \
      ${'--read_structure=' + readStructure} \
      ${'--minimum_quality=' + minimumQuality} \
      ${'--run_start_date=' + runStartDate} \
      $max_reads_in_ram_per_tile \
      $max_records_in_ram \
      --JVMmemory="$mem_in_mb"m \
      $demux_threads \
      ${true='--force_gc=true' false="--force_gc=false" forceGC} \
      --compression_level=5 \
      --loglevel=DEBUG

    illumina.py guess_barcodes --expected_assigned_fraction=0 barcodes.txt metrics.txt barcodes_outliers.txt

    mkdir -p unmatched
    mv Unmatched.bam unmatched/

    OUT_BASENAMES=bam_basenames.txt
    for bam in *.bam; do
      echo "$(basename $bam .bam)" >> $OUT_BASENAMES
    done

    # GNU Parallel refresher:
    # ",," is the replacement string; values after ":::" are substituted where it appears
    parallel --jobs 0 -I ,, \
      "reports.py fastqc \
        ,,.bam \
        ,,_fastqc.html" \
      ::: `cat $OUT_BASENAMES`
  }

  output {
    File        metrics                  = "metrics.txt"
    File        commonBarcodes           = "barcodes.txt"
    File        outlierBarcodes          = "barcodes_outliers.txt"
    Array[File] raw_reads_unaligned_bams = glob("*.bam")
    File        unmatched_reads_bam      = "unmatched/Unmatched.bam"
    Array[File] raw_reads_fastqc         = glob("*_fastqc.html")
    String      viralngs_version         = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "200 GB"
    cpu: 32
    disks: "local-disk 2625 LOCAL"
    dx_instance_type: "mem3_ssd2_v2_x32"
    preemptible: 0  # this is the very first operation before scatter, so let's get it done quickly & reliably
  }
}

task merge_and_reheader_bams {
  Array[File]+  in_bams
  File?         reheader_table # tsv with 3 cols: field, old value, new value
  String        out_basename
  String?       docker="quay.io/broadinstitute/viral-core"

  command {
    set -ex -o pipefail

    read_utils.py --version | tee VERSION

    if [ ${length(in_bams)} -gt 1 ]; then
      read_utils.py merge_bams ${sep=' ' in_bams} merged.bam --loglevel DEBUG
    else
      echo "Skipping merge, only one input file"
      ln -s ${select_first(in_bams)} merged.bam
    fi    

    if [[ -f "${reheader_table}" ]]; then
      read_utils.py reheader_bam merged.bam ${reheader_table} ${out_basename}.bam --loglevel DEBUG
    else
      echo "Skipping reheader, no mapping table specified"
      ln -s merged.bam ${out_basename}.bam
    fi
  }

  output {
    File   out_bam          = "${out_basename}.bam"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "2 GB"
    cpu: 2
    disks: "local-disk 750 LOCAL"
    dx_instance_type: "mem1_ssd2_v2_x4"
  }
}


