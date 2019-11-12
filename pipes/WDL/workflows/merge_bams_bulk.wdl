import "tasks_demux.wdl" as demux

workflow merge_bams_bulk {
    File out_basename_in_bams_table # first column is out_basename, remaining columns are in_bams for that basename
	Array[Array[String]] input_values = read_tsv(out_basename_in_bams_table)

#     Array[File]+ in_bams
#     String out_basename
    String? docker="quay.io/broadinstitute/viral-core"
    
    scatter (input_value in input_values) {
    	call demux.merge_and_reheader_bams {
            input:
            	out_basename = input_value[0],
                in_bams = input_value[1],
                docker = docker
        }
    }
}
