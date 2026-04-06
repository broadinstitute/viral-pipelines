version 1.0

import "../tasks/tasks_metagenomics.wdl" as metagenomics

workflow classify_reads_by_contig {
    meta {
        description: "Classify reads by contig mapping using PAF alignments and VirNucPro contig classifications via DuckDB SQL joins."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File    paf_file
        File    contig_classifications

        Int     min_mapq       = 5
        Float   min_identity   = 90.0
        Float   min_query_cov  = 80.0
    }

    call metagenomics.classify_reads_by_contig as classify_reads {
        input:
            paf_file               = paf_file,
            contig_classifications = contig_classifications,
            min_mapq               = min_mapq,
            min_identity           = min_identity,
            min_query_cov          = min_query_cov
    }

    output {
        File virnucpro_reads_classified = classify_reads.read_classifications
    }
}
