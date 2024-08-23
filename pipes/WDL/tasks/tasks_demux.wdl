version 1.0

task merge_tarballs {
  input {
    Array[File]+ tar_chunks
    String       out_filename

    Int?         machine_mem_gb
    String       docker = "quay.io/broadinstitute/viral-core:2.3.3"
  }

  Int disk_size = 2625

  command {
    set -ex -o pipefail

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi

    file_utils.py --version | tee VERSION

    file_utils.py merge_tarballs \
      ~{out_filename} ~{sep=' ' tar_chunks} \
      --loglevel=DEBUG
  }

  output {
    File   combined_tar     = "~{out_filename}"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu: 16
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd2_v2_x16"
    maxRetries: 2
    preemptible: 0
  }
}

task samplesheet_rename_ids {
  input {
    File   old_sheet
    File?  rename_map
    String old_id_col = 'internal_id'
    String new_id_col = 'external_id'
  }
  parameter_meta { 
    old_sheet: {
      description: "Illumina file with old sample names.",
      category: "required"
    }
    rename_map: {
      description: "New sample name sheet.",
      category: "required"
    }
  }
  String new_base = basename(old_sheet, '.txt')
  Int disk_size = 50
  command <<<
    python3 << CODE
    import csv
    # read in the rename_map file
    old_to_new = {}
    with open('~{default="/dev/null" rename_map}', 'rt') as inf:
      for row in csv.DictReader(inf, delimiter='\t'):
        old_to_new[row['~{old_id_col}']] = row['~{new_id_col}']

    # change all ids in the sample column to new ids
    with open('~{old_sheet}', 'rt') as inf:
      reader = csv.DictReader(inf, delimiter='\t')
      with open('~{new_base}.renamed.txt', 'w', newline='') as outf:
        writer = csv.DictWriter(outf, reader.fieldnames, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for row in reader:
          row['sample'] = old_to_new.get(row['sample'], row['sample'])
          writer.writerow(row)
    CODE
  >>>
  output {
    File new_sheet = '~{new_base}.renamed.txt'
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

task revcomp_i5 {
  input {
    File    old_sheet
    Boolean revcomp=true
    String  docker = "quay.io/broadinstitute/py3-bio:0.1.2"
  }
  String new_base = basename(basename(old_sheet, '.txt'), '.tsv')
  Int disk_size = 50
  command <<<
    python3 << CODE
    import csv
    import Bio.Seq
    old_sheet = "~{old_sheet}"
    revcomp = ~{true="True" false="False" revcomp}

    with open(old_sheet, "rt") as inf:
      with open('~{new_base}'+'.revcompi5.tsv', 'wt') as outf:
        if revcomp:
          sheet = csv.DictReader(inf, delimiter='\t')
          writer = csv.DictWriter(outf, sheet.fieldnames, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
          writer.writeheader()
          for row in sheet:
            if row.get('barcode_2') \
              and all(s in 'ATGC' for s in row['barcode_2']):
              row['barcode_2'] = str(Bio.Seq.Seq(row['barcode_2']).reverse_complement())
            writer.writerow(row)
        else:
          # nop
          outf.write(inf.read())
    CODE
  >>>
  output {
    File new_sheet = '~{new_base}.revcompi5.tsv'
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu:    1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task illumina_demux {
  input {
    File    flowcell_tgz
    Int     lane=1
    Boolean sort_reads=true
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
    Int?    maxRecordsInRam
    Int?    numberOfNegativeControls

    Int?    machine_mem_gb
    Int     disk_size = 2625
    String  docker = "quay.io/broadinstitute/viral-core:2.3.3"
  }

  parameter_meta {
      flowcell_tgz: {
          description: "Illumina BCL directory compressed as tarball. Must contain RunInfo.xml (unless overridden by runinfo), SampleSheet.csv (unless overridden by samplesheet), RTAComplete.txt, and Data/Intensities/BaseCalls/*",
          patterns: ["*.tar.gz", ".tar.zst", ".tar.bz2", ".tar.lz4", ".tgz"]
      }
      samplesheet: {
        description: "CSV file with the library chemistry, sample names and the index tag used for each sample, in addition to some other metrics describing the run.",
        category: "required"
      }
      runinfo: { 
        description: "if we are overriding the RunInfo file, use the path of the file provided. Otherwise the default will be RunInfo.xml. ",
        category: "advanced"
      }
  }

  String out_base = "~{basename(basename(basename(basename(flowcell_tgz, '.zst'), '.gz'), '.tar'), '.tgz')}-L~{lane}"

  command <<<
    set -ex -o pipefail

    # find N% memory
    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 85)

    if [ -z "$TMPDIR" ]; then
      export TMPDIR=$(pwd)
    fi
    FLOWCELL_DIR=$(mktemp -d)

    read_utils.py --version | tee VERSION

    read_utils.py extract_tarball \
      ~{flowcell_tgz} $FLOWCELL_DIR \
      --loglevel=DEBUG

    # if we are overriding the RunInfo file, use the path of the file provided. Otherwise find the file
    if [ -n "~{runinfo}" ]; then
      RUNINFO_FILE="~{runinfo}"
    else
      # full RunInfo.xml path
      RUNINFO_FILE="$(find $FLOWCELL_DIR -type f -name RunInfo.xml | head -n 1)"
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
    if [ "$total_tile_count" -le 2 ]; then
        echo "Detected $total_tile_count tiles, interpreting as MiSeq nano run."
    elif [ "$total_tile_count" -le 8 ]; then
        echo "Detected $total_tile_count tiles, interpreting as MiSeq micro run."
    elif [ "$total_tile_count" -le 16 ]; then
        echo "Detected $total_tile_count tiles, interpreting as iSeq run."
    elif [ "$total_tile_count" -le 28 ]; then
        echo "Detected $total_tile_count tiles, interpreting as MiSeq run."
    elif [ "$total_tile_count" -le 38 ]; then
        echo "Detected $total_tile_count tiles, interpreting as MiSeq run."
    elif [ "$total_tile_count" -le 128 ]; then
        echo "Detected $total_tile_count tiles, interpreting as HiSeq2k run."
    elif [ "$total_tile_count" -le 132 ]; then
        echo "Detected $total_tile_count tiles, interpreting as NextSeq 2000 P2 run."
    elif [ "$total_tile_count" -le 264 ]; then
        echo "Detected $total_tile_count tiles, interpreting as NextSeq 2000 P3 run."
    elif [ "$total_tile_count" -le 288 ]; then
        # increase the number of reads in ram per-tile for NextSeq, since the tiles are larger
        # without this setting, reads will spill to disk and may read the limit
        # on the number of files that can be opened
        # max_reads_in_ram_per_tile=1500000 # deprecared in newer versions of picard, to be removed
        demux_threads=32 # with NovaSeq-size output, OOM errors can sporadically occur with higher thread counts
        max_records_in_ram=1000000
        echo "Detected $total_tile_count tiles, interpreting as NextSeq (mid-output) run."
    elif [ "$total_tile_count" -le 624 ]; then
        demux_threads=32 # with NovaSeq-size output, OOM errors can sporadically occur with higher thread counts
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 80)
        # max_reads_in_ram_per_tile=600000 # deprecared in newer versions of picard, to be removed
        max_records_in_ram=2000000
        echo "Detected $total_tile_count tiles, interpreting as NovaSeq SP run."
    elif [ "$total_tile_count" -le 768 ]; then
        echo "Detected $total_tile_count tiles, interpreting as HiSeq4k run."
    elif [ "$total_tile_count" -le 864 ]; then
        # increase the number of reads in ram per-tile for NextSeq, since the tiles are larger
        # without this setting, reads will spill to disk and may read the limit
        # on the number of files that can be opened
        # max_reads_in_ram_per_tile=1500000 # deprecared in newer versions of picard, to be removed
        max_records_in_ram=2500000
        echo "Detected $total_tile_count tiles, interpreting as NextSeq (high-output) run."
    elif [ "$total_tile_count" -le 896 ]; then
        echo "Detected $total_tile_count tiles, interpreting as HiSeq4k run."
    elif [ "$total_tile_count" -le 1408 ]; then
        demux_threads=32 # with NovaSeq-size output, OOM errors can sporadically occur with higher thread counts
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 80)
        # max_reads_in_ram_per_tile=600000 # deprecared in newer versions of picard, to be removed
        max_records_in_ram=10000000
        echo "Detected $total_tile_count tiles, interpreting as NovaSeq S2 run."
        echo "  **Note: Q20 threshold used since NovaSeq with RTA3 writes only four Q-score values: 2, 12, 23, and 37.**"
        echo "    See: https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/novaseq-hiseq-q30-app-note-770-2017-010.pdf"
    elif [ "$total_tile_count" -le 3744 ]; then
        demux_threads=32 # with NovaSeq-size output, OOM errors can sporadically occur with higher thread counts
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 80)
        max_records_in_ram=2000000
        echo "Detected $total_tile_count tiles, interpreting as NovaSeq run."
        echo "  **Note: Q20 threshold used since NovaSeq with RTA3 writes only four Q-score values: 2, 12, 23, and 37.**"
        echo "    See: https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/novaseq-hiseq-q30-app-note-770-2017-010.pdf"
    elif [ "$total_tile_count" -gt 3744 ]; then
        demux_threads=30 # with NovaSeq-size output, OOM errors can sporadically occur with higher thread counts
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 80)
        # max_reads_in_ram_per_tile=600000 # deprecared in newer versions of picard, to be removed
        max_records_in_ram=2000000
        echo "Tile count: $total_tile_count tiles (unknown instrument type)."
    fi

    # use the passed-in (or default) WDL value first, then fall back to the auto-scaled value
    # if the result of this is null (nothing is passed in, no autoscaled value, no param is passed to the command)
    if [ -n "~{minimumBaseQuality}" ]; then demux_min_base_quality="~{minimumBaseQuality}"; else demux_min_base_quality="$demux_min_base_quality"; fi
    if [ -n "$demux_min_base_quality" ]; then demux_min_base_quality="--minimum_base_quality=$demux_min_base_quality";fi
    
    if [ -n "~{threads}" ]; then demux_threads="~{threads}"; else demux_threads="$demux_threads"; fi
    if [ -n "$demux_threads" ]; then demux_threads="--threads=$demux_threads"; fi
    
    if [ -n "~{maxRecordsInRam}" ]; then max_records_in_ram="~{maxRecordsInRam}"; else max_records_in_ram="$max_records_in_ram"; fi
    if [ -n "$max_records_in_ram" ]; then max_records_in_ram="--max_records_in_ram=$max_records_in_ram"; fi

    # note that we are intentionally setting --threads to about 2x the core
    # count. seems to still provide speed benefit (over 1x) when doing so.
    illumina.py illumina_demux \
      $FLOWCELL_DIR \
      ~{lane} \
      . \
      ~{'--sampleSheet=' + samplesheet} \
      ~{'--runInfo=' + runinfo} \
      ~{'--sequencing_center=' + sequencingCenter} \
      --outMetrics=metrics.txt \
      --commonBarcodes=barcodes.txt \
      ~{'--flowcell=' + flowcell} \
      $demux_min_base_quality \
      ~{'--max_mismatches=' + maxMismatches} \
      ~{'--min_mismatch_delta=' + minMismatchDelta} \
      ~{'--max_no_calls=' + maxNoCalls} \
      ~{'--read_structure=' + readStructure} \
      ~{'--minimum_quality=' + minimumQuality} \
      ~{'--run_start_date=' + runStartDate} \
      ~{true="--sort=true" false="--sort=false" sort_reads} \
      $max_records_in_ram \
      --JVMmemory="$mem_in_mb"m \
      $demux_threads \
      --append_run_id \
      --compression_level=5 \
      --out_meta_by_sample meta_by_sample.json \
      --out_meta_by_filename meta_by_fname.json \
      --out_runinfo runinfo.json \
      --loglevel=DEBUG

    illumina.py guess_barcodes ~{'--number_of_negative_controls ' + numberOfNegativeControls} --expected_assigned_fraction=0 barcodes.txt metrics.txt barcodes_outliers.txt

    illumina.py flowcell_metadata --inDir $FLOWCELL_DIR flowcellMetadataFile.tsv

    mkdir -p unmatched
    mv Unmatched.bam unmatched/

    OUT_BASENAMES=bam_basenames.txt
    for bam in *.bam; do
      echo "$(basename $bam .bam)" >> $OUT_BASENAMES
    done

    # fastqc
    FASTQC_HARDCODED_MEM_PER_THREAD=250 # the value fastqc sets for -Xmx per thread, not adjustable
    num_cpus=$(nproc)
    num_bam_files=$(cat $OUT_BASENAMES | wc -l)
    num_fastqc_jobs=1
    num_fastqc_threads=1
    total_ram_needed_mb=250

    # determine the number of fastqc jobs
    while [[ $total_ram_needed_mb -lt $mem_in_mb ]] && [[ $num_fastqc_jobs -lt $num_cpus ]] && [[ $num_fastqc_jobs -lt $num_bam_files ]]; do
        num_fastqc_jobs=$(($num_fastqc_jobs+1))
        total_ram_needed_mb=$(($total_ram_needed_mb+$FASTQC_HARDCODED_MEM_PER_THREAD))
    done
    # determine the number of fastqc threads per job
    while [[ $(($total_ram_needed_mb)) -lt $mem_in_mb ]] && [[ $(($num_fastqc_jobs*$num_fastqc_threads)) -lt $num_cpus ]]; do
        if [[ $(( $num_fastqc_jobs * $(($num_fastqc_threads+1)) )) -le $num_cpus ]]; then
            num_fastqc_threads=$(($num_fastqc_threads+1))
            total_ram_needed_mb=$(($num_fastqc_jobs*($FASTQC_HARDCODED_MEM_PER_THREAD*$num_fastqc_threads)))
        else
            break
        fi
    done

    # GNU Parallel refresher:
    # ",," is the replacement string; values after ":::" are substituted where it appears
    parallel --jobs $num_fastqc_jobs -I ,, \
      "reports.py fastqc \
        ,,.bam \
        ,,_fastqc.html \
        --out_zip ,,_fastqc.zip \
        --threads $num_fastqc_threads" \
      ::: $(cat $OUT_BASENAMES)

    mv metrics.txt  "~{out_base}-demux_metrics.txt"
    mv runinfo.json "~{out_base}-runinfo.json"

    cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
    cat /proc/loadavg | cut -f 3 -d ' ' > LOAD_15M
    set +o pipefail
    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
  >>>

  output {
    File        metrics                  = "~{out_base}-demux_metrics.txt"
    File        commonBarcodes           = "barcodes.txt"
    File        outlierBarcodes          = "barcodes_outliers.txt"
    Array[File] raw_reads_unaligned_bams = glob("*.bam")
    File        unmatched_reads_bam      = "unmatched/Unmatched.bam"
    Array[File] raw_reads_fastqc         = glob("*_fastqc.html")
    Array[File] raw_reads_fastqc_zip     = glob("*_fastqc.zip")
    Int         max_ram_gb               = ceil(read_float("MEM_BYTES")/1000000000)
    Int         runtime_sec              = ceil(read_float("UPTIME_SEC"))
    Int         cpu_load_15min           = ceil(read_float("LOAD_15M"))

    String      instrument_model         = read_json("~{out_base}-runinfo.json")["sequencer_model"]
    String      flowcell_lane_count      = read_json("~{out_base}-runinfo.json")["lane_count"]

    String      viralngs_version         = read_string("VERSION")

    Map[String,Map[String,String]] meta_by_sample        = read_json('meta_by_sample.json')
    Map[String,Map[String,String]] meta_by_filename      = read_json('meta_by_fname.json')
    Map[String,String]             run_info              = read_json("~{out_base}-runinfo.json")
    File                           meta_by_sample_json   = 'meta_by_sample.json'
    File                           meta_by_filename_json = 'meta_by_fname.json'
    File                           run_info_json         = "~{out_base}-runinfo.json"
  }

  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 200]) + " GB"
    cpu: 32
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem3_ssd2_v2_x32"
    dx_timeout: "20H"
    maxRetries: 2
    preemptible: 0  # this is the very first operation before scatter, so let's get it done quickly & reliably
  }
}

task map_map_setdefault {
  input {
    File          map_map_json
    Array[String] sub_keys
  }
  Int disk_size = 20
  command <<<
    python3 << CODE
    import json
    sub_keys = '~{sep="*" sub_keys}'.split('*')
    with open('~{map_map_json}', 'rt') as inf:
      out = json.load(inf)
    for k in out.keys():
      for sub_key in sub_keys:
        out[k].setdefault(sub_key, "")
    with open('out.json', 'wt') as outf:
      json.dump(out, outf, indent=2)
    CODE
  >>>
  output {
    File out_json = 'out.json'
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task merge_maps {
  input {
    Array[File] maps_jsons
  }
  Int disk_size = 20
  command <<<
    python3 << CODE
    import json
    infiles = '~{sep='*' maps_jsons}'.split('*')
    out = {}
    for fname in infiles:
      with open(fname, 'rt') as inf:
        out.update(json.load(inf))
    with open('out.json', 'wt') as outf:
      json.dump(out, outf, indent=2)
    CODE
  >>>
  output {
    Map[String,Map[String,String]] merged      = read_json('out.json')
    File                           merged_json = 'out.json'
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}
