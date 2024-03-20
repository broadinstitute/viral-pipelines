version 1.0

import "../tasks/tasks_megablast.wdl" as tools

workflow megablast {
    meta {
        desription: "Runs megablast followed by LCA for taxon identification."
        author: "Broad Viral Genomics"
        email: "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }
    input {
        File    inBam
        File    clipDb
        File    blast_db_tgz
        File    taxonomy_db_tgz
        File    taxid_map_file
    }
    call tools.trim_rmdup_subsamp {
        input: 
            inBam = inBam,
            clipDb = clipDb
    }

    call tools.lca_megablast {
        input:
            trimmed_fasta = trim_rmdup_subsamp.trimmed_fasta, 
            blast_db_tgz = blast_db_tgz,
            taxonomy_db_tgz = taxonomy_db_tgz,
            taxdb = taxid_map_file

    }

    output {
        File    LCA_output = lca_megablast.LCA_output
        File    kraken_output_fromat =  lca_megablast.kraken_output_fromat
    }
}