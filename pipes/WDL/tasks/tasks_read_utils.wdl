
task merge_bams {
    Array[File]+    in_bams
    String?         sample_name
    File?           reheader_table # tsv with 3 cols: field, old value, new value
    String          out_basename
    String?         docker="quay.io/broadinstitute/viral-core"

    command {
        set -ex -o pipefail

        read_utils.py --version | tee VERSION
        mem_in_mb=`/opt/viral-ngs/source/docker/calc_mem.py mb 90`

        if [ ${length(in_bams)} -gt 1 ]; then
            read_utils.py merge_bams ${sep=' ' in_bams} merged.bam --JVMmemory="$mem_in_mb"m --loglevel DEBUG
        else
            echo "Skipping merge, only one input file"
            cp ${select_first(in_bams)} merged.bam
        fi    

        if [ -s merged.bam ]; then
            touch reheader_table.txt

            # remap all SM values to user specified value
            if [ -n "${sample_name}" ]; then
              # create sample name remapping table based on existing sample names
              samtools view -H merged.bam | perl -n -e'/SM:(\S+)/ && print "SM\t$1\t'${sample_name}'\n"' | sort | uniq >> reheader_table.txt
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

        else
            # input was empty, so output should be empty (samtools doesn't like empty files)
            touch "${out_basename}.bam"
        fi
    }

    output {
        File   out_bam          = "${out_basename}.bam"
        String viralngs_version = read_string("VERSION")
    }

    runtime {
        docker: "${docker}"
        memory: "3 GB"
        cpu: 4
        disks: "local-disk 750 LOCAL"
        dx_instance_type: "mem1_ssd2_v2_x4"
        preemptible: 0
    }
}

task downsample_bams {
  Array[File]  reads_bam
  Int?         readCount
  Boolean?     deduplicateBefore=false
  Boolean?     deduplicateAfter=false

  String?      docker="quay.io/broadinstitute/viral-core"

  command {
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
        --JVMmemory "1g"
  }

  output {
    Array[File] downsampled_bam  = glob("output/*.downsampled-*.bam")
    String      viralngs_version = read_string("VERSION")
  }
  runtime {
    docker: "${docker}"
    memory: "3 GB"
    cpu:    4
    disks:  "local-disk 750 LOCAL"
    dx_instance_type: "mem1_ssd1_v2_x4"
  }
}


