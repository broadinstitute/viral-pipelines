import "tasks_demux.wdl" as demux

workflow aamerge_bams_bulk {

    Array[File]+ in_bams
    File out_basenames_file # one per line, same order as in_bams_tsv
    String? docker="quay.io/broadinstitute/viral-core"
    
    Array[String] out_basenames = read_lines(out_basenames_file)
    
    scatter (basename_index in range(length(out_basenames))) {
        String out_basename = out_basenames[basename_index]
        
        # identifies the indices of the input bam files containing this output basename
        scatter (in_bams_index in range(length(in_bams))) {
            File in_bam = in_bams[in_bams_index]
            String in_bam_name = basename(in_bam, ".bam")
            
            if(true) {
                File relevant_in_bam = in_bam
            }
        }
        Array[File?] relevant_in_bams = relevant_in_bam # gathers results from the scatter        

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
