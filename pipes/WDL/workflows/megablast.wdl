version 1.0

import "../tasks/tasks_megablast.wdl" as tools

workflow megablast {
    meta {
        desription: "Runs megablast followed by LCA for taxon identification."
        author: "Broad Viral Genomics"
        email: "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    call tools.trim_rmdup_subsamp {
        input: 
            inBam = inBam,
            clipDb = clipDb
    }

    call tools.megablast {
        input:
            cleaned_fasta = trim_rmdup_subsamp.cleaned_fasta, 
            blast_db_tgz = blast_db_tgz,
            taxonomy_db_tgz = taxonomy_db_tgz
    }
    
    output {
        File    LCA_output = megablast.LCA_output
    }
}