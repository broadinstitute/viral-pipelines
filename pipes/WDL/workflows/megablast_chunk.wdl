version 1.0 

import "../tasks/tasks_megablast.wdl" as tools

workflow chunk_megablast {
    meta {
        description: "Chunk megablast function"
        author: "Broad Viral Genomics"
        email: "viral-negs@broadinstitute.org"
        allowNestedInputs: true
    }
    input {
        File    inFasta
        File    blast_db_tgz
        String  db_name 
    }
    call tools.ChunkBlastHits {
        input: 
            inFasta = inFasta, 
            blast_db_tgz = blast_db_tgz,
            db_name = db_name
    }
    output {
        File    blast_hits = ChunkBlastHits.blast_hits
        File    blast_filter_logs = ChunkBlastHits.blast_py_log
        File    chunk_blast_runtime = ChunkBlastHits.duration_seconds
        Int     max_ram_gb = ChunkBlastHits.max_ram_gb
        String  cpu_load = ChunkBlastHits.cpu_load
    }
}
