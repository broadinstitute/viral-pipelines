version 1.0

import "../tasks/task_file_handling.wdl" as file_handling
import "../tasks/task_versioning.wdl" as versioning

workflow bam_to_fastq_pe {
  input {
    File bam_file
    String samplename
  }
  call file_handling.fastq_from_bam_pe {
    input:
    bam_file = bam_file,
    samplename = samplename
  }
  call versioning.version_capture{
    input:
  }
  output {
    String bam_to_fastq_pe_version = version_capture.terra_utilities_version
    String bam_to_fastq_pe_analysis_date = version_capture.date

    File read1 = fastq_from_bam_pe.read1
    File read2 = fastq_from_bam_pe.read2
  }
}
