version 1.0

workflow fastq_to_ubam {
  meta {
    description: "This WDL converts paired FASTQ to uBAM and adds read group information"
    author: "Broad Institute, 2018"
  }
 	call FastqToUBAM
  output {
    File unmapped_bam = FastqToUBAM.unmapped_bam
  }
}

task FastqToUBAM {
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

    Int    machine_mem_gb = 7
    String docker = "broadinstitute/gatk:latest"
  }
  Int command_mem_gb = machine_mem_gb - 1
  Int disk_space_gb = ceil((size(fastq_1, "GB")) * 4 ) + 20
  command {
      /gatk/gatk --java-options "-Xmx~{command_mem_gb}g" \
      FastqToSam \
      --FASTQ ~{fastq_1} \
      ${"--FASTQ2 " + fastq_2}
      --SAMPLE_NAME ~{sample_name} \
      --LIBRARY_NAME ~{library_name} \
      --OUTPUT ~{sample_name}.bam \
      ${"--READ_GROUP_NAME " + readgroup_name} \
      ${"--PLATFORM_UNIT " + platform_unit} \
      ${"--RUN_DATE " + run_date} \
      ${"--PLATFORM " + platform_name} \
      ${"--SEQUENCING_CENTER " + sequencing_center}
  }
  runtime {
    docker: docker
    cpu: 2
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + " HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File unmapped_bam = "~{sample_name}.bam"
  }
}
