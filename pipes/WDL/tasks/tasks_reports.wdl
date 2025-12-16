version 1.0

task alignment_metrics {
  meta {
      description: "Produce various standard metrics and coverage plots via Picard and Samtools for aligned BAM files."
  }

  input {
    File   aligned_bam
    File   ref_fasta
    File?  primers_bed
    String? amplicon_set
    Int?   min_coverage
    Int    max_amp_len=5000
    Int    max_amplicons=500

    Int    machine_mem_gb=16
    String docker = "quay.io/broadinstitute/viral-core:2.5.11"
  }

  String out_basename = basename(aligned_bam, ".bam")
  Int disk_size = 150

  command <<<
    set -e
    MEM_MB=$(free -m | head -2 | tail -1 | awk '{print $4}')
    MEM_MB=$(( MEM_MB > 2048 ? MEM_MB : 2048 ))  # Minimum 2GB heap
    XMX=$(echo "-Xmx"$MEM_MB"m")
    echo "Requesting $MEM_MB MB of RAM for Java"

    # requisite Picard fasta indexing
    python3<<CODE
    import shutil
    import util.file
    with util.file.fastas_with_sanitized_ids("~{ref_fasta}", use_tmp=True) as sanitized_fastas:
        shutil.copyfile(sanitized_fastas[0], 'reference.fasta')
    CODE
    picard $XMX CreateSequenceDictionary -R reference.fasta

    if [ -s "~{ref_fasta}" ]; then
      # get Picard metrics and clean up the junky outputs
      picard $XMX CollectRawWgsMetrics \
        -R reference.fasta \
        -I "~{aligned_bam}" \
        -O picard_raw.raw_wgs_metrics.txt
      grep -v \# picard_raw.raw_wgs_metrics.txt | grep . | head -2 > picard_clean.raw_wgs_metrics.txt

      picard $XMX CollectAlignmentSummaryMetrics \
        -R reference.fasta \
        -I "~{aligned_bam}" \
        -O picard_raw.alignment_metrics.txt
      grep -v \# picard_raw.alignment_metrics.txt | grep . | head -4 > picard_clean.alignment_metrics.txt

      picard $XMX CollectInsertSizeMetrics \
        -I "~{aligned_bam}" \
        -O picard_raw.insert_size_metrics.txt \
        -H picard_raw.insert_size_metrics.pdf \
        --INCLUDE_DUPLICATES true
      grep -v \# picard_raw.insert_size_metrics.txt | grep . | head -2 > picard_clean.insert_size_metrics.txt
    else
      # ref_fasta is empty -> Picard will fail
      touch picard_clean.raw_wgs_metrics.txt picard_clean.alignment_metrics.txt picard_clean.insert_size_metrics.txt
    fi

    # prepend the sample name in order to facilitate tsv joining later
    SAMPLE=$(samtools view -H "~{aligned_bam}" | grep ^@RG | perl -lape 's/^@RG.*SM:(\S+).*$/$1/' | sort | uniq)
    echo -e "sample_sanitized\tbam" > prepend.txt
    echo -e "$SAMPLE\t~{out_basename}" >> prepend.txt
    paste prepend.txt picard_clean.raw_wgs_metrics.txt > "~{out_basename}".raw_wgs_metrics.txt
    echo -e "$SAMPLE\t~{out_basename}" >> prepend.txt
    echo -e "$SAMPLE\t~{out_basename}" >> prepend.txt
    paste prepend.txt picard_clean.alignment_metrics.txt > "~{out_basename}".alignment_metrics.txt
    echo -e "sample_sanitized\tbam" > prepend.txt
    echo -e "$SAMPLE\t~{out_basename}" >> prepend.txt
    paste prepend.txt picard_clean.insert_size_metrics.txt > "~{out_basename}".insert_size_metrics.txt

    touch "~{out_basename}".ampliconstats.txt "~{out_basename}".ampliconstats_parsed.txt
    echo -e "sample_sanitized\tbam\tamplicon_set\tamplicon_idx\tamplicon_left\tamplicon_right\tFREADS\tFDEPTH\tFPCOV\tFAMP" > "~{out_basename}.ampliconstats_parsed.txt"
    if [ -n "~{primers_bed}" ]; then
      # samtools ampliconstats
      cat "~{primers_bed}" | sort -k 1,1 -k 4,4 -t $'\t' > primers-sorted_for_samtools.bed
      set +e # there are just some weird bed files out there -- let them fail silently
      samtools ampliconstats -s -@ $(nproc) \
        ~{'-d ' + min_coverage} \
        ~{'-l ' + max_amp_len} \
        ~{'-a ' + max_amplicons} \
        -o "~{out_basename}".ampliconstats.txt primers-sorted_for_samtools.bed "~{aligned_bam}"

      # parse into our own tsv to facilitate tsv joining later
      if [ -n "~{default='' amplicon_set}" ]; then
        AMPLICON_SET="~{default='' amplicon_set}"
      else
        AMPLICON_SET=$(basename "~{primers_bed}" .bed)
      fi
      grep ^AMPLICON "~{out_basename}".ampliconstats.txt | cut -f 2- > AMPLICON
      grep ^FREADS "~{out_basename}".ampliconstats.txt | cut -f 3- | tr '\t' '\n' > FREADS; echo "" >> FREADS
      grep ^FDEPTH "~{out_basename}".ampliconstats.txt | cut -f 3- | tr '\t' '\n' > FDEPTH; echo "" >> FDEPTH
      grep ^FPCOV  "~{out_basename}".ampliconstats.txt | cut -f 3- | tr '\t' '\n' > FPCOV;  echo "" >> FPCOV
      grep ^FAMP   "~{out_basename}".ampliconstats.txt | cut -f 4 | tail +2 > FAMP
      for i in $(cut -f 1 AMPLICON); do echo -e "$SAMPLE\t~{out_basename}\t$AMPLICON_SET"; done > prepend.txt
      wc -l prepend.txt AMPLICON FREADS FDEPTH FPCOV FAMP
      paste prepend.txt AMPLICON FREADS FDEPTH FPCOV FAMP | grep '\S' >> "~{out_basename}.ampliconstats_parsed.txt"
    fi
  >>>

  output {
    File wgs_metrics         = "~{out_basename}.raw_wgs_metrics.txt"
    File alignment_metrics   = "~{out_basename}.alignment_metrics.txt"
    File insert_size_metrics = "~{out_basename}.insert_size_metrics.txt"
    File amplicon_stats      = "~{out_basename}.ampliconstats.txt"
    File amplicon_stats_parsed = "~{out_basename}.ampliconstats_parsed.txt"
  }

  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem3_ssd1_v2_x4"
    maxRetries: 2
  }
}

task plot_coverage {
  input {
    File    aligned_reads_bam
    String  sample_name

    Boolean skip_mark_dupes           = false
    Boolean plot_only_non_duplicates  = false
    Boolean bin_large_plots           = false
    String? binning_summary_statistic = "max" # max or min

    Int? plot_width_pixels            = 1100
    Int? plot_height_pixels           = 850
    Int? plot_pixels_per_inch         = 100

    Int? max_coverage_depth
    Int? base_q_threshold
    Int? mapping_q_threshold
    Int? read_length_threshold
    String? plotXLimits # of the form "min max" (ints, space between)
    String? plotYLimits # of the form "min max" (ints, space between)

    String  docker = "quay.io/broadinstitute/viral-core:2.5.11"
  }

  Int disk_size = 375
  
  command {
    set -ex -o pipefail

    read_utils.py --version | tee VERSION

    samtools view -c ${aligned_reads_bam} | tee reads_aligned
    if [ "$(cat reads_aligned)" != "0" ]; then
      samtools index -@ "$(nproc)" "${aligned_reads_bam}"

      PLOT_DUPE_OPTION=""
      if [[ "${skip_mark_dupes}" != "true" ]]; then
        PLOT_DUPE_OPTION="${true='--plotOnlyNonDuplicates' false="" plot_only_non_duplicates}"
      fi
      
      BINNING_OPTION="${true='--binLargePlots' false="" bin_large_plots}"

      # plot coverage
      reports.py plot_coverage \
        "${aligned_reads_bam}" \
        "${sample_name}.coverage_plot.pdf" \
        --outSummary "${sample_name}.coverage_plot.txt" \
        --plotFormat pdf \
        --plotWidth 1100 \
        --plotHeight 850 \
        --plotDPI 100 \
        ${"-m " + max_coverage_depth} \
        ${"-q " + base_q_threshold} \
        ${"-Q " + mapping_q_threshold} \
        ${"-l " + read_length_threshold} \
        ${"--plotXLimits " + plotXLimits} \
        ${"--plotYLimits " + plotYLimits} \
        $PLOT_DUPE_OPTION \
        $BINNING_OPTION \
        --binningSummaryStatistic ${binning_summary_statistic} \
        --plotTitle "${sample_name} coverage plot" \
        --loglevel=DEBUG

    else
      touch ${sample_name}.coverage_plot.pdf ${sample_name}.coverage_plot.txt
    fi

    # collect figures of merit
    set +o pipefail # grep will exit 1 if it fails to find the pattern
    samtools view -H ${aligned_reads_bam} | perl -n -e'/^@SQ.*LN:(\d+)/ && print "$1\n"' |  python -c "import sys; print(sum(int(x) for x in sys.stdin))" | tee assembly_length
    # report only primary alignments 260=exclude unaligned reads and secondary mappings
    samtools view -h -F 260 ${aligned_reads_bam} | samtools flagstat - | tee ${sample_name}.flagstat.txt
    grep properly ${sample_name}.flagstat.txt | cut -f 1 -d ' ' | tee read_pairs_aligned
    samtools view ${aligned_reads_bam} | cut -f10 | tr -d '\n' | wc -c | tee bases_aligned
    python -c "print (float("$(cat bases_aligned)")/"$(cat assembly_length)") if "$(cat assembly_length)">0 else print(0)" > mean_coverage
  }

  output {
    File   coverage_plot      = "${sample_name}.coverage_plot.pdf"
    File   coverage_tsv       = "${sample_name}.coverage_plot.txt"
    Int    assembly_length    = read_int("assembly_length")
    Int    reads_aligned      = read_int("reads_aligned")
    Int    read_pairs_aligned = read_int("read_pairs_aligned")
    Float  bases_aligned      = read_float("bases_aligned")
    Float  mean_coverage      = read_float("mean_coverage")
    String viralngs_version   = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "7 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x4"
    preemptible: 1
    maxRetries: 2
  }
}

task merge_coverage_per_position {
  input {
    Array[File]+ coverage_tsvs
    File         ref_fasta

    String       out_report_name = "coverage_report.csv"
    Int          disk_size = 100
    String       docker = "quay.io/broadinstitute/py3-bio:0.1.2"
  }

  command <<<
    set -e

    python3<<CODE
    import os
    import pandas as pd
    from functools import reduce
    import Bio.SeqIO

    # get genome length
    genome_length = 0
    with open('~{ref_fasta}', 'rt') as inf:
      for seq in Bio.SeqIO.parse(inf, 'fasta'):
        genome_length += len(seq.seq.replace("-",""))

    # Loop through a list of file paths and read in each depth.tsv generated as part of assemble_refbased
    depths_dfs = []
    for in_tsv in ("~{sep='", "' coverage_tsvs}"):
        sample_name = '.'.join(os.path.basename(in_tsv).split('.')[:-2])
        sample_depths_df = pd.read_csv(in_tsv, sep='\t', header=None
            ).rename(columns={0:'Ref',1:'Position',2:sample_name})
        depths_dfs.append(sample_depths_df)

    # Condense all depths into a single dataframe
    df_merged = reduce(lambda left,right:
        pd.merge(left,right,on=['Ref','Position'],how='outer'),
        depths_dfs)
    df_merged = df_merged.fillna(0)

    #Create dummy df that contains all positions along the genome
    dummy_df = pd.DataFrame([range(1,genome_length)]).T.rename(columns={0:'Position'})
    df_merged = df_merged.merge(dummy_df, on='Position', how='right').fillna(0)
    df_merged = df_merged.drop(['Ref'], axis=1)
    df_merged.to_csv("~{out_report_name}", index=False)
    CODE
  >>>

  output {
    File   coverage_multi_sample_per_position_csv  = out_report_name
  }

  runtime {
    docker: "${docker}"
    memory: "2 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd2_v2_x4"
    maxRetries: 2
  }
}

task coverage_report {
  input {
    Array[File]+ mapped_bams
    Array[File]  mapped_bam_idx = []  # optional.. speeds it up if you provide it, otherwise we auto-index
    String       out_report_name = "coverage_report.txt"

    String       docker = "quay.io/broadinstitute/viral-core:2.5.11"
  }

  Int disk_size = 375

  command <<<
    set -e
    reports.py --version | tee VERSION
    python3 << CODE
    import tools.samtools
    import reports
    samtools = tools.samtools.SamtoolsTool()
    in_bams = list([bam for bam in ["~{sep='", "' mapped_bams}"] if bam and not samtools.isEmpty(bam)])
    if in_bams:
      reports.coverage_only(in_bams, "~{out_report_name}")
    else:
      with open("~{out_report_name}", "w") as outf:
        outf.write('\t'.join(('sample', 'aln2self_cov_median', 'aln2self_cov_mean', 'aln2self_cov_mean_non0', 'aln2self_cov_1X', 'aln2self_cov_5X', 'aln2self_cov_20X', 'aln2self_cov_100X'))+'\n')
    CODE
  >>>

  output {
    File   coverage_report  = out_report_name
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd2_v2_x4"
    maxRetries: 2
  }
}

task assembly_bases {
    meta {
      description: "Count bases in a fasta file."
    }

    input {
      File   fasta
      String docker ="ubuntu"
    }

    Int disk_size = 50

    command {
        set -e
        grep -v '^>' "~{fasta}" | tr -d '\n' | wc -c | tee assembly_length
        grep -v '^>' "~{fasta}" | tr -d '\nNn' | wc -c | tee assembly_length_unambiguous
    }

    output {
        Int assembly_length             = read_int("assembly_length")
        Int assembly_length_unambiguous = read_int("assembly_length_unambiguous")
    }

    runtime {
        docker: "${docker}"
        memory: "1 GB"
        cpu: 1
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
}

task fastqc {
  input {
    File   reads_bam

    String docker = "quay.io/broadinstitute/viral-core:2.5.11"
  }
  parameter_meta {
    reads_bam:{ 
    description: "Input reads in BAM format.",
    category: "required"
    }

  }

  String   reads_basename=basename(reads_bam, ".bam")
  Int disk_size = 375

  command {
    set -ex -o pipefail
    reports.py --version | tee VERSION
    reports.py fastqc ${reads_bam} ${reads_basename}_fastqc.html --out_zip ${reads_basename}_fastqc.zip
  }

  output {
    File   fastqc_html      = "${reads_basename}_fastqc.html"
    File   fastqc_zip       = "${reads_basename}_fastqc.zip"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: "2 GB"
    cpu: 1
    docker: "${docker}"
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task align_and_count {
  input {
    File   reads_bam
    File   ref_db
    Int    topNHits = 3

    Boolean filter_bam_to_proper_primary_mapped_reads         = true
    Boolean do_not_require_proper_mapped_pairs_when_filtering = false
    Boolean keep_singletons_when_filtering                    = false
    Boolean keep_duplicates_when_filtering                    = false

    Int?   cpu
    Int?   machine_mem_gb
    String docker = "quay.io/broadinstitute/viral-core:2.5.11"
  }

  String  reads_basename=basename(reads_bam, ".bam")
  String  ref_basename=basename(ref_db, ".fasta")
  Int disk_size = ceil((10 * size(reads_bam, "GB") + 2 * size(ref_db, "GB") + 150) / 375.0) * 375

  # Autoscale CPU based on input size: 4 CPUs for small inputs, up to 16 CPUs for larger inputs
  # Linear scaling: 4 + (input_GB / 15) * 60, capped at 16, rounded to nearest multiple of 4
  # NOTE: Capped low because minimap2_idxstats doesn't parallelize well - see broadinstitute/viral-core#145
  Float        cpu_unclamped = 4.0 + (size(reads_bam, "GB") / 15.0) * 60.0
  Int          cpu_actual = select_first([cpu, floor(((if cpu_unclamped > 16.0 then 16.0 else cpu_unclamped) + 2.0) / 4.0) * 4])
  # Memory scales with CPU at 2x ratio (default), or use override
  Int          machine_mem_gb_actual = select_first([machine_mem_gb, cpu_actual * 2])

  parameter_meta {
    reads_bam: {
      description: "Unaligned reads in BAM format",
      pattern: ["*.bam"],
      category: "required"
    }
    ref_db: {
      description: "Reference genome in FASTA format",
      pattern: ["*.FASTA"],
      category: "required"
    }
    filter_bam_to_proper_primary_mapped_reads: {
      description: "If specified, reads till be filtered after alignment to include only those flagged as properly paired.",
      category: "optional"
    }
    do_not_require_proper_mapped_pairs_when_filtering: {
      description: "Do not require reads to be properly paired when filtering",
      category: "optional"
    }
    keep_singletons_when_filtering: {
      description: "Keep singletons when filtering",
      category: "optional"
    }
    keep_duplicates_when_filtering: {
      description: "Do not exclude reads marked as duplicates when filtering",
      category: "optional"
    }
  }
  command <<<
    set -ex -o pipefail

    read_utils.py --version | tee VERSION

    ln -s "~{reads_bam}" "~{reads_basename}.bam"
    read_utils.py minimap2_idxstats \
      "~{reads_basename}.bam" \
      "~{ref_db}" \
      --outStats "~{reads_basename}.count.~{ref_basename}.txt.unsorted" \
      ~{true="--filterReadsAfterAlignment"   false="" filter_bam_to_proper_primary_mapped_reads} \
      ~{true="--doNotRequirePairsToBeProper" false="" do_not_require_proper_mapped_pairs_when_filtering} \
      ~{true="--keepSingletons"              false="" keep_singletons_when_filtering} \
      ~{true="--keepDuplicates"              false="" keep_duplicates_when_filtering} \
      --loglevel=DEBUG

    sort -b -r -n -k3 "~{reads_basename}.count.~{ref_basename}.txt.unsorted" > "~{reads_basename}.count.~{ref_basename}.txt"
    head -n ~{topNHits} "~{reads_basename}.count.~{ref_basename}.txt" > "~{reads_basename}.count.~{ref_basename}.top_~{topNHits}_hits.txt"
    TOP_HIT="$(head -1 '~{reads_basename}.count.~{ref_basename}.txt' | cut -f 1 | sed 's/\*/\\*/' | tee '~{reads_basename}.count.~{ref_basename}.top.txt')"

    TOTAL_COUNT_OF_TOP_HIT=$(grep -E "^($TOP_HIT)" "~{reads_basename}.count.~{ref_basename}.txt" | cut -f3 | tee TOTAL_COUNT_OF_TOP_HIT)
    TOTAL_COUNT_OF_LESSER_HITS=$((grep -vE "^(\*|$TOP_HIT)" "~{reads_basename}.count.~{ref_basename}.txt" || echo "0" ) | cut -f3 | paste -sd+ - | bc -l | tee TOTAL_COUNT_OF_LESSER_HITS)
    echo $TOTAL_COUNT_OF_TOP_HIT | tee TOTAL_COUNT_OF_TOP_HIT
    echo $TOTAL_COUNT_OF_LESSER_HITS | tee TOTAL_COUNT_OF_LESSER_HITS

    if [ $TOTAL_COUNT_OF_LESSER_HITS -ne 0 -o $TOTAL_COUNT_OF_TOP_HIT -ne 0 ]; then
      PCT_MAPPING_TO_LESSER_HITS=$( echo "scale=3; 100 * $TOTAL_COUNT_OF_LESSER_HITS / ($TOTAL_COUNT_OF_LESSER_HITS + $TOTAL_COUNT_OF_TOP_HIT)" | \
        bc -l | awk '{printf "%.3f\n", $0}' | tee '~{reads_basename}.count.~{ref_basename}.pct_lesser_hits_of_mapped.txt' )
    else
      echo "PCT_MAPPING_TO_LESSER_HITS cannot be calculated: there were no hits to any sequence"
      PCT_MAPPING_TO_LESSER_HITS=$( echo "null" | tee '~{reads_basename}.count.~{ref_basename}.pct_lesser_hits_of_mapped.txt')
    fi

    TOTAL_READS_IN_INPUT=$(samtools view -c "~{reads_basename}.bam")
    echo $TOTAL_READS_IN_INPUT | tee TOTAL_READS_IN_INPUT
    if [ $TOTAL_READS_IN_INPUT -eq 0 ]; then
      echo "no reads in input bam"
      PCT_OF_INPUT_READS_MAPPED=$(echo "0" | tee "~{reads_basename}.count.~{ref_basename}.pct_total_reads_mapped.txt")
      echo "PCT_TOP_HIT_OF_TOTAL_READS cannot be calculated: there were no mapping hits, or no reads"
      PCT_TOP_HIT_OF_TOTAL_READS=$( echo "null" | tee '~{reads_basename}.count.~{ref_basename}.pct_top_hit_of_total_reads.txt')
    else
      PCT_OF_INPUT_READS_MAPPED=$( echo "scale=3; 100 * ($TOTAL_COUNT_OF_LESSER_HITS + $TOTAL_COUNT_OF_TOP_HIT) / $TOTAL_READS_IN_INPUT" | \
      bc -l | awk '{printf "%.3f\n", $0}' | tee '~{reads_basename}.count.~{ref_basename}.pct_total_reads_mapped.txt' )

      if [ $TOTAL_COUNT_OF_TOP_HIT -ne 0 ]; then
        PCT_TOP_HIT_OF_TOTAL_READS=$( echo "scale=3; 100 * ($TOTAL_COUNT_OF_TOP_HIT / $TOTAL_READS_IN_INPUT)" | \
          bc -l | awk '{printf "%.3f\n", $0}' | tee '~{reads_basename}.count.~{ref_basename}.pct_top_hit_of_total_reads.txt' )
      else
        echo "PCT_TOP_HIT_OF_TOTAL_READS cannot be calculated: there were no mapping hits, or no reads"
        PCT_TOP_HIT_OF_TOTAL_READS=$( echo "null" | tee '~{reads_basename}.count.~{ref_basename}.pct_top_hit_of_total_reads.txt')
      fi
    fi

    cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
    cat /proc/loadavg | cut -f 3 -d ' ' > LOAD_15M
    set +o pipefail
    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi; } > MEM_BYTES
  >>>

  output {
    File   report           = "~{reads_basename}.count.~{ref_basename}.txt"
    
    File   report_top_hits  = "~{reads_basename}.count.~{ref_basename}.top_~{topNHits}_hits.txt"
    String top_hit_id       = read_string("~{reads_basename}.count.~{ref_basename}.top.txt")

    Int    reads_total          = read_int("TOTAL_READS_IN_INPUT")
    Int    reads_mapped_top_hit = read_int("TOTAL_COUNT_OF_TOP_HIT")
    Int    reads_mapped         = read_int("TOTAL_COUNT_OF_LESSER_HITS") + read_int("TOTAL_COUNT_OF_TOP_HIT")

    String pct_total_reads_mapped     = read_string('~{reads_basename}.count.~{ref_basename}.pct_total_reads_mapped.txt')
    String pct_top_hit_of_total_reads = read_string('~{reads_basename}.count.~{ref_basename}.pct_top_hit_of_total_reads.txt')
    String pct_lesser_hits_of_mapped  = read_string('~{reads_basename}.count.~{ref_basename}.pct_lesser_hits_of_mapped.txt')

    Int    max_ram_gb       = ceil(read_float("MEM_BYTES")/1000000000)
    Int    runtime_sec      = ceil(read_float("UPTIME_SEC"))
    Int    cpu_load_15min   = ceil(read_float("LOAD_15M"))
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: machine_mem_gb_actual + " GB"
    cpu: cpu_actual
    docker: "${docker}"
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x4"
    maxRetries: 2
  }
}

task align_and_count_summary {
  input {
    Array[File]+ counts_txt

    String       output_prefix = "count_summary"

    String       docker = "quay.io/broadinstitute/viral-core:2.5.11"
  }

  Int disk_size = 100

  command {
    set -ex -o pipefail

    reports.py --version | tee VERSION
    reports.py aggregate_alignment_counts ${sep=' ' counts_txt} "${output_prefix}".tsv --loglevel=DEBUG
  }

  output {
    File   count_summary    = "${output_prefix}.tsv"
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: "7 GB"
    cpu: 8
    docker: "${docker}"
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task aggregate_metagenomics_reports {
  input {
    Array[File]+ kraken_summary_reports 
    String       aggregate_taxon_heading_space_separated  = "Viruses"
    String       aggregate_taxlevel_focus                 = "species"
    Int          aggregate_top_N_hits                     = 5

    String       docker = "quay.io/broadinstitute/viral-classify:2.5.1.0"
  }

  parameter_meta {
    aggregate_taxon_heading_space_separated: { description: "The taxonomic heading to analyze. More than one can be specified." }
    aggregate_taxlevel_focus:                { description: "species,genus,family,order,class,phylum,kingdom,superkingdom" }
    aggregate_top_N_hits:                    { description: "only include the top N hits from a given sample in the aggregate report" }
  }

  String       aggregate_taxon_heading = sub(aggregate_taxon_heading_space_separated, " ", "_") # replace spaces with underscores for use in filename
  Int disk_size = 50

  command {
    set -ex -o pipefail

    metagenomics.py --version | tee VERSION
    metagenomics.py taxlevel_summary \
      ${sep=' ' kraken_summary_reports} \
      --csvOut aggregate_taxa_summary_${aggregate_taxon_heading}_by_${aggregate_taxlevel_focus}_top_${aggregate_top_N_hits}_by_sample.csv \
      --noHist \
      --taxHeading ${aggregate_taxon_heading_space_separated} \
      --taxlevelFocus ${aggregate_taxlevel_focus} \
      --zeroFill --includeRoot --topN ${aggregate_top_N_hits} \
      --loglevel=DEBUG
  }

  output {
    File   krakenuniq_aggregate_taxlevel_summary = "aggregate_taxa_summary_${aggregate_taxon_heading}_by_${aggregate_taxlevel_focus}_top_${aggregate_top_N_hits}_by_sample.csv"
    String viralngs_version                      = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: "3 GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd2_v2_x2"
    preemptible: 0
    maxRetries: 2
  }
}

task MultiQC {
  input {
    Array[File]    input_files

    String?        title
    String?        comment
    String?        file_name
    String         out_dir       = "./multiqc-output"
    String?        template
    String?        tag
    String?        ignore_analysis_files
    String?        ignore_sample_names
    File?          sample_names
    Array[String]? exclude_modules
    Array[String]? module_to_use
    String?        output_data_format
    Boolean        force         = false
    Boolean        full_names    = false
    Boolean        data_dir      = false
    Boolean        no_data_dir   = false
    Boolean        zip_data_dir  = false
    Boolean        export        = false
    Boolean        flat          = false
    Boolean        interactive   = true
    Boolean        lint          = false
    Boolean        pdf           = false
    Boolean        megaQC_upload = false # Upload generated report to MegaQC if MegaQC options are found
    File?          config  # directory
    String?        config_yaml

    String         docker = "quay.io/biocontainers/multiqc:1.32--pyhdfd78af_1"
  }

  parameter_meta {
    output_data_format: { description: "[tsv|yaml|json] default:tsv" }
  }

  # get the basename in all wdl use the filename specified (sans ".html" extension, if specified)
  String report_filename = if (defined(file_name)) then basename(select_first([file_name]), ".html") else "multiqc"
  Int disk_size = 375

  command {
      set -ex -o pipefail

      echo "${sep='\n' input_files}" > input-filenames.txt
      echo "" >> input-filenames.txt

      # Run MultiQC but allow it to fail (it crashes with exit 1 on empty/invalid zip files)
      set +e
      multiqc \
      --file-list input-filenames.txt \
      --outdir "${out_dir}" \
      ${true="--force" false="" force} \
      ${true="--fullnames" false="" full_names} \
      ${"--title " + title} \
      ${"--comment " + comment} \
      ${"--filename " + file_name} \
      ${"--template " + template} \
      ${"--tag " + tag} \
      ${"--ignore " + ignore_analysis_files} \
      ${"--ignore-samples" + ignore_sample_names} \
      ${"--sample-names " + sample_names} \
      ${true="--exclude " false="" defined(exclude_modules)}${sep=' --exclude ' select_first([exclude_modules,[]])} \
      ${true="--module " false="" defined(module_to_use)}${sep=' --module ' select_first([module_to_use,[]])} \
      ${true="--data-dir" false="" data_dir} \
      ${true="--no-data-dir" false="" no_data_dir} \
      ${"--data-format " + output_data_format} \
      ${true="--zip-data-dir" false="" zip_data_dir} \
      ${true="--export" false="" export} \
      ${true="--flat" false="" flat} \
      ${true="--interactive" false="" interactive} \
      ${true="--lint" false="" lint} \
      ${true="--pdf" false="" pdf} \
      ${false="--no-megaqc-upload" true="" megaQC_upload} \
      ${"--config " + config} \
      ${"--cl-config " + config_yaml }
      MULTIQC_EXIT_CODE=$?
      set -e

      # Ensure output directory exists (MultiQC may remove it if no results found or if it crashed)
      mkdir -p "${out_dir}"

      if [ -z "${file_name}" ]; then
        mv "${out_dir}/${report_filename}_report.html" "${out_dir}/${report_filename}.html" 2>/dev/null || true
      fi

      # Create placeholder HTML report if MultiQC didn't create one (happens when no valid results found or on crash)
      if [ ! -f "${out_dir}/${report_filename}.html" ]; then
        echo "<!DOCTYPE html><html><head><meta charset=\"UTF-8\"><title>MultiQC Report</title></head><body><h1>MultiQC Report</h1><p>No analysis results found in input files.</p></body></html>" > "${out_dir}/${report_filename}.html"
      fi

      # Ensure data directory exists before tarring (MultiQC only creates it when results are found)
      mkdir -p "${out_dir}/${report_filename}_data"
      tar -c "${out_dir}/${report_filename}_data" | gzip -c > "${report_filename}_data.tar.gz"
  }

  output {
      File multiqc_report           = "${out_dir}/${report_filename}.html"
      File multiqc_data_dir_tarball = "${report_filename}_data.tar.gz"
  }

  runtime {
    memory: "8 GB"
    cpu: 16
    docker: "${docker}"
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem2_ssd1_v2_x2"
    maxRetries: 2
  }
}

task compare_two_genomes {
  input {
    File   genome_one
    File   genome_two
    String out_basename

    String docker = "quay.io/broadinstitute/viral-assemble:2.5.1.0"
  }

  Int disk_size = 50

  command <<<
    set -ex -o pipefail
    assembly.py --version | tee VERSION
    assembly.py alignment_summary "~{genome_one}" "~{genome_two}" --outfileName "~{out_basename}.txt" --printCounts --loglevel=DEBUG
    cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
    cat /proc/loadavg > CPU_LOAD
    set +o pipefail
    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
  >>>

  output {
    File   comparison_table = "~{out_basename}.txt"
    Int    max_ram_gb       = ceil(read_float("MEM_BYTES")/1000000000)
    Int    runtime_sec      = ceil(read_float("UPTIME_SEC"))
    String cpu_load         = read_string("CPU_LOAD")
    String viralngs_version = read_string("VERSION")
  }

  runtime {
    memory: "3 GB"
    cpu: 2
    docker: docker
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 1
    maxRetries: 2
  }
}


