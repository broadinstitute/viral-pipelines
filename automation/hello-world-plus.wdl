version 1.0

workflow HelloWorldPlus {

	input {
    	File input_file
    	String id
    }
  
  	call echoFileName { 
  		input:
        	input_file = input_file,
            output_base_name = id
	}
   
    output {
    	File output_file = echoFileName.output_file
    }
}

task echoFileName {

	input {
		File input_file
    	String output_base_name
    }
    
    String base_file_name=basename(input_file)
    String output_file_name="~{output_base_name}_output.txt"

	command <<<
        echo "The name of the input file was ~{base_file_name}." > ~{output_file_name}
    >>>

    runtime {
    	docker: "ubuntu"
        cpu: 1
        disks: "local-disk 1 HDD"
        preemptible: 1
    }
    output {
    	File output_file = output_file_name
    }
}