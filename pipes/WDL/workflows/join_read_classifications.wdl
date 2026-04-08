version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow join_read_classifications {
    meta {
        description: "Join read-level classifications from Kallisto, Kraken2, VirNucPro, and geNomad into a single ZSTD-compressed Parquet file keyed on SAMPLE_ID + READ_ID."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File?   kallisto_summary
        File?   kraken2_reads
        File?   vnp_reads
        File?   genomad_virus_summary
        File?   centrifuger_reads
        String  sample_id
        Boolean filter_human_only_k2 = true
    }

    call metagenomics.join_read_classifications as join_reads {
        input:
            kallisto_summary      = kallisto_summary,
            kraken2_reads         = kraken2_reads,
            vnp_reads             = vnp_reads,
            genomad_virus_summary = genomad_virus_summary,
            centrifuger_reads     = centrifuger_reads,
            sample_id             = sample_id,
            filter_human_only_k2  = filter_human_only_k2
    }

    output {
        File classifications_summary = join_reads.classifications_parquet
    }
}
