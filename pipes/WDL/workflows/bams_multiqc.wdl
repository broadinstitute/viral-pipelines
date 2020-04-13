import "../tasks/tasks_reports.wdl" as reports

workflow bams_multiqc {
    Array[File]+  read_bams

    scatter(reads_bam in read_bams) {
        call reports.fastqc as fastqc {
            input:
                reads_bam = reads_bam
        }
    }

    call reports.MultiQC {
        input:
            input_files = fastqc.fastqc_zip
    }

}
