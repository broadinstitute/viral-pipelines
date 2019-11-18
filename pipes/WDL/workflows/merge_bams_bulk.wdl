import "tasks_demux.wdl" as demux

workflow merge_bams_bulk {
    Map[String, Array[File]] map_of_strings_to_arrays_of_files
    Array[Pair] array_of_pairs
    Pair[String, File] string_and_file_pair
    Array[Pair[File, String]] array_of_pairs_of_files_and_strings
    Map[File, String] map_of_files_to_strings
    Array[Array[File]] array_of_arrays_of_files
    Array[Pair[String, Array[File]]] array_of_pairs_of_strings_and_arrays_of_files
    Array[Array[File]] array_of_arrays_of_files
    Array[String] array_of_strings
    
    File? reheader_table
    String? docker="quay.io/broadinstitute/viral-core"
}

#     Array[File]+ in_bams
#     File out_basenames # one per line
#     File? reheader_table
#     String? docker="quay.io/broadinstitute/viral-core"
#     
#     # identifies out_basename for each in_bam file
#     Array[String] out_basenames_list = read_lines(out_basenames)
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
