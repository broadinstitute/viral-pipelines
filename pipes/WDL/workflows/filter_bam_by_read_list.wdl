version 1.0

import "../tasks/tasks_read_utils.wdl" as read_utils

workflow filter_bam_by_read_list {
    meta {
        description: "Subset a BAM to only reads listed in a flat READ_ID text file, with mate-level specificity. Each line is a READ_ID optionally ending in '/1' or '/2'; bare READ_IDs keep both mates."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File    reads_bam
        File    read_ids_txt
        String? sample_id_override
    }

    call read_utils.filter_bam_by_read_list as filter_bam {
        input:
            reads_bam          = reads_bam,
            read_ids_txt       = read_ids_txt,
            sample_id_override = sample_id_override
    }

    output {
        File filtered_bam = filter_bam.filtered_bam
    }
}
