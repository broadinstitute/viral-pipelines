import "tasks_demux.wdl" as demux

workflow merge_bams_bulk {
    Array[File]+ in_bams # any order
    File in_bam_out_bam_table # first column: input bam file basename, second column: output bam file basename (one line per INPUT file)
    File? reheader_table
    String? docker="quay.io/broadinstitute/viral-core"
    
    # removes ".bam"s from ends of filenames in in_bam_out_bam_table
    call clean_in_bam_out_bam_table {
        input: table = in_bam_out_bam_table
    }
    
    # generates map with key: input bam file name -> value: output bam file basename
    Map[String, String] in_bam_to_out_bam = read_map(clean_in_bam_out_bam_table.clean_table)
    
    # retrieves unique output bam file basenames (no repeats)
    call unique_values_in_second_column {
        input: table = clean_in_bam_out_bam_table.clean_table
    }
    
    # collects and merges input bam files for each output bam file
    scatter (out_bam in unique_values_in_second_column.unique_values) {
        String placeholder = write_map(in_bam_to_out_bam) # need to touch map in outer scatter for it to be seen in inner scatter
        
        # retrieves the input bam files for this output bam file
        scatter (in_bam in in_bams) {
            String in_bam_basename = basename(in_bam, ".bam")
            if(in_bam_to_out_bam[in_bam_basename] == out_bam) {
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

task clean_in_bam_out_bam_table {
    File table
    
    command {
    	cat ${table} | sed 's/[.]bam$//g' | sed $'s/[.]bam\t/\t/g' | tee in_bam_out_bam_table
    }
    
    output {
        File clean_table = "in_bam_out_bam_table"
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
