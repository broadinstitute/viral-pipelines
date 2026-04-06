version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow summarize_kb_extract_reads {
    meta {
        description: "Summarize kb extract read output with taxonomy annotation. Extracts reads from tarball, joins with taxonomy mapping, and outputs zstd-compressed TSV."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File    extract_reads_tar
        File    taxonomy_map_csv
        String  taxonomy_level = "highest"
    }

    call metagenomics.summarize_kb_extract_reads as summarize {
        input:
            extract_reads_tar = extract_reads_tar,
            taxonomy_map_csv = taxonomy_map_csv,
            taxonomy_level = taxonomy_level
    }

    output {
        File kallisto_read_classifications = summarize.summary_tsv_zst
    }
}
