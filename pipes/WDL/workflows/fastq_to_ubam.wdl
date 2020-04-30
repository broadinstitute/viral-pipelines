version 1.0

# TO DO: convert to single task with conditional logic embedded in task command, utilize viral-core read_utils.py fastq_to_bam

workflow ConvertFastQsToUnmappedBamWf {
  meta {
    description: "This WDL converts paired FASTQ to uBAM and adds read group information"
    author: "Broad Institute, 2018"
  }

  input {
    File   fastq_1 
    File?  fastq_2 

    String sample_name 
    String readgroup_name 
    String library_name 
    String platform_unit 
    String run_date 
    String platform_name 
    String sequencing_center 
  }

  # Convert pair of FASTQs to uBAM
  if(defined(fastq_2)) {
  	call PairedFastQsToUnmappedBAM {
      input:
        sample_name = sample_name,
        fastq_1 = fastq_1,
        fastq_2 = fastq_2,
        readgroup_name = readgroup_name,
        library_name = library_name,
        platform_unit = platform_unit,
        run_date = run_date,
        platform_name = platform_name,
        sequencing_center = sequencing_center
    }
  }
  
  if(!defined(fastq_2)) {
  	call SingleFastQsToUnmappedBAM {
      input:
        sample_name = sample_name,
        fastq_1 = fastq_1,
        readgroup_name = readgroup_name,
        library_name = library_name,
        platform_unit = platform_unit,
        run_date = run_date,
        platform_name = platform_name,
        sequencing_center = sequencing_center
    }
  }
  
  # Outputs that will be retained when execution is complete
  output {
    File? paired_unmapped_bam = PairedFastQsToUnmappedBAM.paired_output_unmapped_bam
    File? single_unmapped_bam = SingleFastQsToUnmappedBAM.single_output_unmapped_bam
  }
}


task PairedFastQsToUnmappedBAM {
  input {
    # Command parameters
    String sample_name
    File fastq_1
    File? fastq_2
    String readgroup_name
    String library_name
    String platform_unit
    String run_date
    String platform_name
    String sequencing_center

    # Runtime parameters
    Int addtional_disk_space_gb = 10
    Int machine_mem_gb = 7
    Int preemptible_attempts = 3
    String docker = "broadinstitute/gatk:latest"
  }
    Int command_mem_gb = machine_mem_gb - 1
    Int disk_space_gb = ceil((size(fastq_1, "GB") + size(fastq_2, "GB")) * 2 ) + addtional_disk_space_gb
  
  command {
      /gatk/gatk --java-options "-Xmx~{command_mem_gb}g" \
      FastqToSam \
      --FASTQ ~{fastq_1} \
      --FASTQ2 ~{fastq_2} \
      --OUTPUT ~{sample_name}.bam \
      --READ_GROUP_NAME ~{readgroup_name} \
      --SAMPLE_NAME ~{sample_name} \
      --LIBRARY_NAME ~{library_name} \
      --PLATFORM_UNIT ~{platform_unit} \
      --RUN_DATE ~{run_date} \
      --PLATFORM ~{platform_name} \
      --SEQUENCING_CENTER ~{sequencing_center}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
  }
  output {
    File paired_output_unmapped_bam = "~{sample_name}.bam"
  }
}

task SingleFastQsToUnmappedBAM {
  input {
    # Command parameters
    String sample_name
    File fastq_1
    String readgroup_name
    String library_name
    String platform_unit
    String run_date
    String platform_name
    String sequencing_center

    # Runtime parameters
    Int addtional_disk_space_gb = 10
    Int machine_mem_gb = 7
    Int preemptible_attempts = 3
    String docker = "broadinstitute/gatk:latest"
  }
    Int command_mem_gb = machine_mem_gb - 1
    Int disk_space_gb = ceil((size(fastq_1, "GB")) * 2 ) + addtional_disk_space_gb
  
  command {
      /gatk/gatk --java-options "-Xmx~{command_mem_gb}g" \
        FastqToSam \
        --FASTQ ~{fastq_1} \
        --OUTPUT ~{sample_name}.bam \
        --READ_GROUP_NAME ~{readgroup_name} \
        --SAMPLE_NAME ~{sample_name} \
        --LIBRARY_NAME ~{library_name} \
        --PLATFORM_UNIT ~{platform_unit} \
        --RUN_DATE ~{run_date} \
        --PLATFORM ~{platform_name} \
        --SEQUENCING_CENTER ~{sequencing_center}
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
  }
  output {
    File single_output_unmapped_bam = "~{sample_name}.bam"
  }
}
