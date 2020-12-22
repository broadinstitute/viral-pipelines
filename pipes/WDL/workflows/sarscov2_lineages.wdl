version 1.0

import "../tasks/tasks_sarscov2.wdl" as sarscov2

workflow sarscov2_lineages {
    meta {
        description: "Call Nextclade and Pangolin lineages on a single SARS-CoV-2 genome"
    }

    input {
    	File genome_fasta
    }

    call sarscov2.nextclade_one_sample {
    	input:
    		genome_fasta = genome_fasta
    }

    call sarscov2.pangolin_one_sample {
    	input:
    		genome_fasta = genome_fasta
    }

    output {
    	String nextclade_clade = nextclade_one_sample.nextclade_clade
    	File   nextclade_tsv   = nextclade_one_sample.nextclade_tsv
    	String pangolin_clade  = pangolin_one_sample.pangolin_clade
    	File   pangolin_csv    = pangolin_one_sample.pangolin_csv
    }
}
