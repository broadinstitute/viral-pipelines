import "tasks_demux.wdl" as demux

workflow aamerge_bams_bulk {

    Array[File]+ in_bams
    File out_basenames # one per line
    File? reheader_table
    String? docker="quay.io/broadinstitute/viral-core"
    
    Array[String] out_basenames_list = read_lines(out_basenames)
    scatter (out_basename in out_basenames_list) {
        
        # identifies the indices of the input bam files containing this output basename
        scatter (in_bam in in_bams) {
            call does_in_bam_match_out_basename {
                input:
                    out_basename = out_basename,
                    in_bam = in_bam
            }
            
            if(does_in_bam_match_out_basename.match) {
                File relevant_in_bam = in_bam
            }
        }
        Array[File?] relevant_in_bams_optional = relevant_in_bam # gathers results from the scatter        Array[File] relevant_in_bams = select_all(relevant_in_bams_optional)
        Array[File] relevant_in_bams = select_all(relevant_in_bams_optional)

        # merges the bam files to produce this output file
        call demux.merge_and_reheader_bams {
            input:
                out_basename = out_basename,
                in_bams = relevant_in_bams,
                reheader_table = reheader_table,
                docker = docker
        }
    }
}

# returns true if the basename of in_bam contains out_basename,
# either at the start or end of the string or surrounded by any of [._-]
task does_in_bam_match_out_basename {
    String out_basename
    File in_bam
    
    String in_bam_name = basename(in_bam, ".bam")
    
    command {
        # basename (exact match)
        if [[ ${in_bam_name} =~ ^${out_basename}$ ]]; then
            echo true | tee match
        # something[._-]basename
        elif [[ ${in_bam_name} =~ [._-]${out_basename}$ ]]; then
            echo true | tee match
        # basename[._-]something
        elif [[ ${in_bam_name} =~ ^${out_basename}[._-] ]]; then
            echo true | tee match
        # something[._-]basename[._-]something
        elif [[ ${in_bam_name} =~ [._-]${out_basename}[._-] ]]; then
            echo true | tee match
        else
            echo false | tee match
        fi
    }
    
    output {
        Boolean match = read_boolean("match")
    }
}