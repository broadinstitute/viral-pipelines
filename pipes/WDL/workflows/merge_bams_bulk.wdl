import "tasks_demux.wdl" as demux

workflow merge_bams_bulk {
    Array[File]+ in_bams # any order
    File in_bam_out_bam_table # first column: input bam file basename, second column: output bam file basename
    File out_bams_file
    File? reheader_table
    String? docker="quay.io/broadinstitute/viral-core"
    
    # generates map with key: input bam file name -> value: output bam file basename
    Map[String, String] in_bam_to_out_bam = read_map(in_bam_out_bam_table)
    
    # retrieves unique output bam file basenames (no repeats)
#     call unique_values_in_second_column {
#         input: table = in_bam_out_bam_table
#     }
#     Array[String] out_bams = unique_values_in_second_column.unique_values
    Array[String] out_bams = read_lines(out_bams_file)
    
    # collects and merges input bam files for each output bam file
    scatter (out_bam in out_bams) {
        # retrieves the input bam files for this output bam file
        scatter (in_bam in in_bams) {
            String in_bam_basename_long = basename(in_bam)
            String in_bam_basename_short = basename(in_bam, ".bam")
            if(defined(in_bam_to_out_bam[in_bam_basename_long]) && in_bam_to_out_bam[in_bam_basename_long] == out_bam
              || defined(in_bam_to_out_bam[in_bam_basename_short]) && in_bam_to_out_bam[in_bam_basename_short] == out_bam) {
                File relevant_in_bam = in_bam
            }
        }
        Array[File] relevant_in_bams = select_all(relevant_in_bam) # gathers input bam files from the scatter
        
        # merges the relevant input bam files to produce this output file
        call demux.merge_and_reheader_bams {
            input:
                out_basename = basename(out_bam, ".bam"),
                in_bams = relevant_in_bams,
                reheader_table = reheader_table,
                docker = docker
        }
    }
}

task unique_values_in_second_column {
    File table
    
    command {
        cut -f2 ${table} | sort | uniq | tee unique_values
    }
    
    output {
        Array[String] unique_values = read_lines("unique_values")
    }
}

# workflow merge_bams_bulk {
# 
#     Array[File]+ in_bams
#     File out_basenames # one per line
#     File? reheader_table
#     String? docker="quay.io/broadinstitute/viral-core"
#     
#     # identifies out_basename for each in_bam file
#     Array[String] out_basenames_list = read_lines(out_basenames)
#     String test = out_basenames_list[0]
# #     Map[String, Int] test = {"a": 1, "b": 2}
# #     Int fun = test["b"]
# #     Map[File, String] in_bam_to_out_basename = {}
# #     scatter (in_bam in in_bams) {
# #         call get_out_basename {
# #             input:
# #                 in_bam = in_bam,
# #                 out_basenames = out_basenames_list
# #         }
# #         String out_basename = get_out_basename.out_basename
# #         in_bam_to_out_basename[in_bam] = out_basename
# #     }
#     
#     # generates an output file for each out_basename
#     scatter (out_basename in out_basenames_list) {
#         
#         # identifies the input bam files containing this output basename
#         # (surrounded by start or end of string or any of [._-])
#         scatter (in_bam in in_bams) {
#             if(in_bam_to_out_basename[in_bam] == out_basename) {
#                 File relevant_in_bam = in_bam
#             }
#         }
#         Array[File] relevant_in_bams = select_all(relevant_in_bam) # gathers results from the scatter
# 
#         # merges the relevant input bam files to produce this output file
#         call demux.merge_and_reheader_bams {
#             input:
#                 out_basename = out_basename,
#                 in_bams = relevant_in_bams,
#                 reheader_table = reheader_table,
#                 docker = docker
#         }
#     }
# }
# 
# task get_out_basename {
#     File in_bam
#     Array[String] out_basenames
#     
#     String in_bam_basename = basename(in_bam, ".bam")
#     
#     command {
#         for out_basename in ${sep=' ' out_basenames}; do
#             # basename (exact match)
#             if [[ ${in_bam_name} =~ ^$out_basename$ ]]; then
#                 echo true | tee out_basename
#             # something[._-]basename
#             elif [[ ${in_bam_name} =~ [._-]$out_basename$ ]]; then
#                 echo true | tee out_basename
#             # basename[._-]something
#             elif [[ ${in_bam_name} =~ ^$out_basename[._-] ]]; then
#                 echo true | tee out_basename
#             # something[._-]basename[._-]something
#             elif [[ ${in_bam_name} =~ [._-]$out_basename[._-] ]]; then
#                 echo true | tee out_basename
#             else
#                 echo false | tee out_basename
#             fi
#         done
#     }
#     output {
#         String out_basename = read_string("out_basename")
#     }
# }

#     scatter (out_basename in read_lines(out_basenames)) {
#         
#         # identifies the input bam files containing this output basename
#         # (surrounded by start or end of string or any of [._-])
#         scatter (in_bam in in_bams) {
#             call does_in_bam_match_out_basename {
#                 input:
#                     out_basename = out_basename,
#                     in_bam = in_bam
#             }
#             
#             if(does_in_bam_match_out_basename.match) {
#                 File relevant_in_bam = in_bam
#             }
#         }
#         Array[File] relevant_in_bams = select_all(relevant_in_bam) # gathers results from the scatter
# 
#         # merges the relevant input bam files to produce this output file
#         call demux.merge_and_reheader_bams {
#             input:
#                 out_basename = out_basename,
#                 in_bams = relevant_in_bams,
#                 reheader_table = reheader_table,
#                 docker = docker
#         }
#     }
# 
# # returns true if the basename of in_bam contains out_basename,
# # separated from other characters by start or end of string or any of [._-]
# task does_in_bam_match_out_basename {
#     String out_basename
#     File in_bam
#     
#     String in_bam_name = basename(in_bam, ".bam")
#     
#     command {
#         # basename (exact match)
#         if [[ ${in_bam_name} =~ ^${out_basename}$ ]]; then
#             echo true | tee match
#         # something[._-]basename
#         elif [[ ${in_bam_name} =~ [._-]${out_basename}$ ]]; then
#             echo true | tee match
#         # basename[._-]something
#         elif [[ ${in_bam_name} =~ ^${out_basename}[._-] ]]; then
#             echo true | tee match
#         # something[._-]basename[._-]something
#         elif [[ ${in_bam_name} =~ [._-]${out_basename}[._-] ]]; then
#             echo true | tee match
#         else
#             echo false | tee match
#         fi
#     }
#     
#     output {
#         Boolean match = read_boolean("match")
#     }
# }
