version 1.0

import "../tasks/tasks_megablast.wdl" as tools
import "../tasks/tasks_metagenomics.wdl" as metagenomics

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
        File    krona_taxonomy_tab
        File    taxonomy_db_tgz
        Int     host_species
        String  db_name
        String sample_name = basename(in_bam, '.bam')
    }
    parameter_meta {
        krona_taxonomy_tab: {
            description: "Krona taxonomy database containing a single file: 'taxonomy.tab' (exact name), or possibly just a compressed 'taxonomy.tab"
            patterns: ["*.tab.zst", "*.tab.gz", "*.tab", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
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
    call metagenomics.krona {
        input: 
            reports_txt_gz = blastoff.blastoff_kraken,
            krona_taxonomy_db_tgz = krona_taxonomy_tab,
            input_type = "tsv",
            out_basename = "~{sample_name}.krona"
    }

    }
    output {
        File    most_popular_taxon_id = blastoff.most_popular_taxon_id
        File    blastoff_txt_results = blastoff.blastoff_results
        File    blastoff_kraken = blastoff.blastoff_kraken
        File    krona_html = metagenomics.krona_report_html
    }
}