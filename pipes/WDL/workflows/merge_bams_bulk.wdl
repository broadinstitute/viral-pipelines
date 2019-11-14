import "tasks_demux.wdl" as demux

workflow merge_bams_bulk {

    Array[File]+ in_bams
    File out_basenames_file # one per line, same order as in_bams_tsv
    String? docker="quay.io/broadinstitute/viral-core"
	
    Array[String] out_basenames = read_lines(out_basenames_file)
    
    scatter (basename_index in range(length(out_basenames))) {
        # retrieves output file basename and list of input filenames for this row
        String out_basename = out_basenames[basename_index]
        
        # identifies the indices of the input bam files containing this output basename
#         scatter (in_bams_index in range(length(in_bams))) {
#             in_bam = in_bams[in_bams_index]
#             in_bam_name = basename(in_bam, ".bam")
#             
#             if(true) {
#                 relevant_in_bam_index = in_bams_index
#             }
#         }
#         Array[Int] relevant_in_bam_indices = relevant_in_bam_index # gathers results from the scatter 

        Array[Int] relevant_in_bam_indices = range(length(in_bams))
        
        # retrieves the input bam files corresponding to the filenames
#         scatter (relevant_in_bam_index in relevant_in_bam_indices) {
#             relevant_in_bam = in_bams[relevant_in_bam_index]
#         }
#         Array[File] relevant_in_bams = relevant_in_bam # gathers results from the scatter
        
        Array[File] relevant_in_bams = in_bams
        
        # merges the bam files to produce this output file
        call demux.merge_and_reheader_bams {
            input:
            	out_basename = out_basename,
                in_bams = relevant_in_bams,
                docker = docker
        }
    }
}

# File in_bams_tsv # filenames to merge, tab-separated, each line one output file
# 	Array[Array[File]] in_bams_table = read_tsv(in_bams_tsv)

# 	Array[Array[String]] in_bams_filenames_table = read_tsv(in_bams_tsv)
#         Array[String] relevant_in_bams_filenames = in_bams_filenames_table[basename_index]

# goes through an array of files and identifies the indices of those that contain a pattern
# task subset_files_list {
#     Array[File]+ files
#     String pattern
#     
#     scatter (index in range(length(out_basenames))) {
#     }
#     
#     output {
#     }
# }

        # retrieves the input bam files corresponding to the filenames
#         scatter(relevant_in_bams_filename in relevant_in_bams_filenames) {
#             call subset_files_list {
#                 input:
#                     input_files = in_bams,
#                     exact_pattern = relevant_in_bams_filename
#             }
#             relevant_in_bam = subset_files_list.matching_files
#         }
#         Array[File] relevant_in_bams = relevant_in_bam # gathers results from the scatter

# task subset_files_list {
#     Array[File]+ input_files
#     String pattern
#     
#     output {
#         Array[File] all_files = input_files
#         Array[File] matching_files = glob(pattern)
#     }
# }

        
#         scatter(in_bam in in_bams) {
#             in_bam_filename = basename(in_bam, ".bam")
#             
#             # input bam files matching the input table
#             # TODO
#             
#             # input bam files containing this output basename
#             # TODO
#             
#             relevant_in_bam = in_bam
#         }
#         Array[File] relevant_in_bams = relevant_in_bam # gathers results from the scatter



#     File out_basename_in_bams_table # first column is out_basename, remaining columns are in_bams for that basename
# 	Array[Array[String]] input_values = read_tsv(out_basename_in_bams_table)
# 
#     scatter (input_value in input_values) {
#         String out_basename = input_value[0]
#         Array[String] in_bam_filenames = [input_value[1], input_value[2]]
#     	call demux.merge_and_reheader_bams {
#             input:
#             	out_basename = out_basename,
#                 in_bams = in_bams,
#                 docker = docker
#         }
#     }
# 
# 
# 
#     File out_basenames_file
#     Array[String] out_basenames = read_lines(out_basenames_file)
#     
#     Array[File]+ in_bams
#     String? docker="quay.io/broadinstitute/viral-core"
#     
#     scatter (out_basename in out_basenames) {
#         these_in_bams = in_bams
#         
#         if()
#     }

#     String   sample_name = basename(basename(basename(reads_unmapped_bam, ".bam"), ".taxfilt"), ".clean")

# }

# task subset_files_list_subset {
#   Array[File]+ files
#   Array[String]+ included_filenames
# 
#   command {
#     set -ex -o pipefail
# 
#     read_utils.py --version | tee VERSION
# 
#     if [ ${length(in_bams)} -gt 1 ]; then
#       read_utils.py merge_bams ${sep=' ' in_bams} merged.bam --loglevel DEBUG
#     else
#       echo "Skipping merge, only one input file"
#       ln -s ${select_first(in_bams)} merged.bam
#     fi    
# 
#     if [[ -f "${reheader_table}" ]]; then
#       read_utils.py reheader_bam merged.bam ${reheader_table} ${out_basename}.bam --loglevel DEBUG
#     else
#       echo "Skipping reheader, no mapping table specified"
#       ln -s merged.bam ${out_basename}.bam
#     fi
#   }
# 
#   output {
#     File   out_bam          = "${out_basename}.bam"
#   }
# 
#   runtime {
#     docker: "${docker}"
#     memory: "2000 MB"
#     cpu: 2
#     dx_instance_type: "mem1_ssd2_v2_x4"
#   }
# }