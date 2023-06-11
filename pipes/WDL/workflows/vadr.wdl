version 1.0 

import "../tasks/tasks_vadr.wdl" as vadr_kit

workflow vadr {
    meta {

        description: "VADR test run"
        author: "fnegrete"
        email: "viral_ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    call vadr_kit.vadr_tool{
        input: 
            genome_fasta = genome_fasta,
            library_accession_no = library_accession_no
    }
}