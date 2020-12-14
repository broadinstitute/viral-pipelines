version 1.0

task merge_and_reheader_bams {
    meta {
      description: "Merge and/or reheader bam files using a mapping table. This task can modify read group tags in a BAM header field for single BAM files or as part of a BAM merge operation. The output is a single BAM file (given one or more input BAMs) and a three-column tab delimited text table that defines: the field, the old value, and the new value (e.g. LB, old_lib_name, new_lib_name or SM, old_sample_name, new_sample_name)"
    }

    input {
      Array[File]+    in_bams
      String?         sample_name
      File?           reheader_table
      String          out_basename

      String          docker="quay.io/broadinstitute/viral-core:2.1.12"
    }

    command {
        set -ex -o pipefail

        read_utils.py --version | tee VERSION
        mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

        if [ ${length(in_bams)} -gt 1 ]; then
            read_utils.py merge_bams ${sep=' ' in_bams} merged.bam --JVMmemory="$mem_in_mb"m --loglevel DEBUG
        else
            echo "Skipping merge, only one input file"
            cp ${sep=' ' in_bams} merged.bam
        fi    

        # remap all SM values to user specified value
        if [ -n "${sample_name}" ]; then
          # create sample name remapping table based on existing sample names
          samtools view -H merged.bam | perl -n -e'/SM:(\S+)/ && print "SM\t$1\t'"${sample_name}"'\n"' | sort | uniq >> reheader_table.txt
        fi

        # remap arbitrary headers using user specified table
        if [[ -f "${reheader_table}" ]]; then
          cat "${reheader_table}" >> reheader_table.txt
        fi

        # reheader bam file if requested
        if [ -s reheader_table.txt ]; then
          read_utils.py reheader_bam merged.bam reheader_table.txt "${out_basename}.bam" --loglevel DEBUG
        else
          mv merged.bam "${out_basename}.bam"
        fi
    }

    output {
        File   out_bam          = "${out_basename}.bam"
        String viralngs_version = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: "3 GB"
        cpu: 2
        disks: "local-disk 750 LOCAL"
        dx_instance_type: "mem1_ssd2_v2_x4"
        preemptible: 0
    }
}

task rmdup_ubam {
  meta {
    description: "Perform read deduplication on unaligned reads."
  }

  input {
    File     reads_unmapped_bam
    String   method="mvicuna"

    Int?     machine_mem_gb
    String?  docker="quay.io/broadinstitute/viral-core:2.1.12"
  }

  parameter_meta {
    reads_unmapped_bam: { description: "unaligned reads in BAM format", patterns: ["*.bam"] }
    method:             { description: "mvicuna or cdhit" }
  }

  String reads_basename = basename(reads_unmapped_bam, ".bam")
  
  command {
    set -ex -o pipefail
    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)
    read_utils.py --version | tee VERSION

    read_utils.py rmdup_"${method}"_bam \
      "${reads_unmapped_bam}" \
      "${reads_basename}".dedup.bam \
      --JVMmemory "$mem_in_mb"m \
      --loglevel=DEBUG

    samtools view -c ${reads_basename}.dedup.bam | tee dedup_read_count_post
    reports.py fastqc ${reads_basename}.dedup.bam ${reads_basename}.dedup_fastqc.html --out_zip ${reads_basename}.dedup_fastqc.zip
  }

  output {
    File   dedup_bam             = "${reads_basename}.dedup.bam"
    File   dedup_fastqc          = "${reads_basename}.dedup_fastqc.html"
    File   dedup_fastqc_zip      = "${reads_basename}.dedup_fastqc.zip"
    Int    dedup_read_count_post = read_int("dedup_read_count_post")
    String viralngs_version      = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 7]) + " GB"
    cpu:    2
    disks:  "local-disk 375 LOCAL"
    dx_instance_type: "mem2_ssd1_v2_x2"
  }
}

task downsample_bams {
  meta {
    description: "Downsample reads in a BAM file randomly subsampling to a target read count. Read deduplication can occur either before or after random subsampling, or not at all (default: not at all)."
  }

  input {
    Array[File]+  reads_bam
    Int?         readCount
    Boolean?     deduplicateBefore=false
    Boolean?     deduplicateAfter=false

    Int?         machine_mem_gb
    String       docker="quay.io/broadinstitute/viral-core:2.1.12"
  }

  command {
    set -ex -o pipefail

    # find 90% memory
    mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

    if [[ "${deduplicateBefore}" == "true" ]]; then
      DEDUP_OPTION="--deduplicateBefore"
    elif [[ "${deduplicateAfter}" == "true" ]]; then
      DEDUP_OPTION="--deduplicateAfter"
    fi

    if [[ "${deduplicateBefore}" == "true" && "${deduplicateAfter}" == "true" ]]; then
      echo "deduplicateBefore and deduplicateAfter are mutually exclusive. Only one can be used."
      exit 1
    fi
    
    read_utils.py --version | tee VERSION

    read_utils.py downsample_bams \
        ${sep=' ' reads_bam} \
        --outPath ./output \
        ${'--readCount=' + readCount} \
        $DEDUP_OPTION \
        --JVMmemory "$mem_in_mb"m
  }

  output {
    Array[File] downsampled_bam  = glob("output/*.downsampled-*.bam")
    String      viralngs_version = read_string("VERSION")
  }

  runtime {
    docker: "${docker}"
    memory: select_first([machine_mem_gb, 3]) + " GB"
    cpu:    4
    disks:  "local-disk 750 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x4"
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
    String? platform_name
    String? sequencing_center

    String  docker="quay.io/broadinstitute/viral-core:2.1.12"
  }
  parameter_meta {
    fastq_1: { description: "Unaligned read1 file in fastq format", patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"] }
    fastq_2: { description: "Unaligned read2 file in fastq format. This should be empty for single-end read conversion and required for paired-end reads. If provided, it must match fastq_1 in length and order.", patterns: ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"] }
    sample_name: { description: "Sample name. This is required and will populate the 'SM' read group value and will be used as the output filename (must be filename-friendly)." }
    library_name: { description: "Library name. This is required and will populate the 'LB' read group value. SM & LB combinations must be identical for any sequencing reads generated from the same sequencing library, and must be distinct for any reads generated from different libraries." }
  }
  command {
      set -ex -o pipefail

      # find 90% memory
      mem_in_mb=$(/opt/viral-ngs/source/docker/calc_mem.py mb 90)

      read_utils.py --version | tee VERSION

      picard -Xmx"$mem_in_mb"m \
        FastqToSam \
        FASTQ="~{fastq_1}" \
        ${"FASTQ2=" + fastq_2} \
        SAMPLE_NAME="${sample_name}" \
        LIBRARY_NAME="${library_name}" \
        OUTPUT="${sample_name}".bam \
        ${"READ_GROUP_NAME=" + readgroup_name} \
        ${"PLATFORM_UNIT=" + platform_unit} \
        ${"RUN_DATE=" + run_date} \
        ${"PLATFORM=" + platform_name} \
        ${"SEQUENCING_CENTER=" + sequencing_center}
  }
  runtime {
    docker: docker
    cpu: 2
    memory: "3 GB"
    disks: "local-disk 375 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File unmapped_bam = "~{sample_name}.bam"
  }
}
