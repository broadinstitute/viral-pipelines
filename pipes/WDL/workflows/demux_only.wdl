import "../tasks/tasks_demux.wdl" as tasks_demux
import "../tasks/tasks_reports.wdl" as reports

workflow demux_only {
    call tasks_demux.illumina_demux

    call reports.MultiQC {
        input:
            input_files = illumina_demux.raw_reads_fastqc
    }
}
