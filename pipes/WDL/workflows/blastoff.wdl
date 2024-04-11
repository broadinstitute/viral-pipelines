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
        Int     db_name
    }
    call tools.trim_rmdup_subsamp {
        input: 
            inBam = inBam,
            clipDb = clipDb
    }
    call tools.blastoff {
        input:
            trimmed_fasta = trim_rmdup_subsamp.trimmed_fasta, 
            blast_db_tgz = blast_db_tgz,
            taxonomy_db_tgz = taxonomy_db_tgz,
            host_species = host_species,
            db_name = db_name

    }
    output {
        File    most_popular_taxon_id = blastoff.most_popular_taxon_id
    }
}