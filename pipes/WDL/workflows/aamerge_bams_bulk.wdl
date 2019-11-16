import "tasks_demux.wdl" as demux

workflow aamerge_bams_bulk {

    Array[File]+ in_bams
    File out_basenames # one per line
    File? in_bam_basenames # tab-separated; one line per output file listed in out_basenames
    File? reheader_table
    String? docker="quay.io/broadinstitute/viral-core"
    
    Array[String] out_basenames_list = read_lines(out_basenames)
    scatter (out_basename in out_basenames_list) {
        
        if(defined(in_bam_basenames)) {
            # we have an input file listing the input bams for each output bam,
            # so we will use it
            Array[File] relevant_in_bams = in_bams
        }
        
        if(!defined(in_bam_basenames) == false) {
            # there is no input file listing the input bams for each output bam,
            # so we will merge all potential input bams matching the output bam
            # basename at start or end or surrounded by any of [._-]
            
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
            Array[File?] relevant_in_bams_optional = relevant_in_bam # gathers results from the scatter
            Array[File] relevant_in_bams = select_all(relevant_in_bams_optional)
        }

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
    File in_bam
    String out_basename
    
    String in_bam_name = basename(in_bam, ".bam")

    command {
        if [[ ${in_bam_name} =~ ^${out_basename}$ ]]; then
            echo true | tee match
        elif [[ ${in_bam_name} =~ [._-]${out_basename}$ ]]; then
            echo true | tee match
        elif [[ ${in_bam_name} =~ ^${out_basename}[._-] ]]; then
            echo true | tee match
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

# returns true if the basename of in_bam exactly matches expected_in_bam_name
task does_in_bam_match_expected_in_bam {
    File in_bam
    String expected_in_bam_name
    
    String in_bam_name = basename(in_bam, ".bam")
    
    command {
        if [[ ${in_bam_name} =~ ^${expected_in_bam_name}$ ]]; then
            echo true | tee match
        else
            echo false | tee match
        fi
    }

    output {
        Boolean match = read_boolean("match")
    }
}

#     Array[Int] basename_scatter_range = range(length(out_basenames_list))
#     scatter (basename_index in basename_scatter_range) {
#         String out_basename = out_basenames_list[basename_index]

# File in_bams_tsv # filenames to merge, tab-separated, each line one output file
#     Array[Array[File]] in_bams_table = read_tsv(in_bams_tsv)

#     Array[Array[String]] in_bams_filenames_table = read_tsv(in_bams_tsv)
#         Array[String] relevant_in_bams_filenames = in_bams_filenames_table[basename_index]

# task subset_files_list {
#     Array[File]+ input_files
#     String pattern
#     
#     output {
#         Array[File] all_files = input_files
#         Array[File] matching_files = glob(pattern)
#     }
# }

#     String   sample_name = basename(basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt"), ".clean")
