import "tasks_demux.wdl" as demux

workflow merge_bams_bulk {
    Array[File]+ in_bams
    File? reheader_table # tsv with 3 cols: field, old value, new value
    String out_basename
    String? docker="quay.io/broadinstitute/viral-core"

    call demux.merge_and_reheader_bams {
        input:
            in_bams = in_bams,
            reheader_table = reheader_table,
            out_basename = out_basename,
            docker = docker
    }
}
