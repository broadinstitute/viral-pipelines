version 1.0

import "../tasks/tasks_reports.wdl" as reports

workflow diff_genome_sets {

    input {
        Array[File]   genome_set_one
        Array[File]   genome_set_two
    }

    scatter(sample in zip(genome_set_one, genome_set_two)) {
        call reports.compare_two_genomes {
            input:
                genome_one = sample.left,
                genome_two = sample.right,
                out_basename = basename(sample.left, '.fasta')
        }
    }

    call reports.tsv_stack {
        input:
            input_tsvs = compare_two_genomes.comparison_table,
            out_basename = "diff_genome_sets"
    }

    output {
        File diff = tsv_stack.out_tsv
    }

}
