version 1.1

import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_sarscov2.wdl" as sarscov2

workflow sarscov2_lineages {
    meta {
        description: "Call Nextclade and Pangolin lineages on a single SARS-CoV-2 genome"
    }

    input {
        File genome_fasta
    }

    call nextstrain.nextclade_one_sample {
        input:
            genome_fasta = genome_fasta,
            dataset_name  = "sars-cov-2"
    }

    call sarscov2.pangolin_one_sample {
        input:
            genome_fasta = genome_fasta
    }

    output {
        String nextclade_clade    = nextclade_one_sample.nextclade_clade
        File   nextclade_tsv      = nextclade_one_sample.nextclade_tsv
        File   nextclade_json     = nextclade_one_sample.nextclade_json
        String nextclade_aa_subs  = nextclade_one_sample.aa_subs_csv
        String nextclade_aa_dels  = nextclade_one_sample.aa_dels_csv
        String nextclade_version  = nextclade_one_sample.nextclade_version
        String pango_lineage      = pangolin_one_sample.pango_lineage
        String pangolin_conflicts = pangolin_one_sample.pangolin_conflicts
        String pangolin_notes     = pangolin_one_sample.pangolin_notes
        String scorpio_call       = pangolin_one_sample.scorpio_call
        String scorpio_notes      = pangolin_one_sample.scorpio_notes
        File   pango_lineage_report = pangolin_one_sample.pango_lineage_report
        String pangolin_docker    = pangolin_one_sample.pangolin_docker
        String pangolin_versions  = pangolin_one_sample.pangolin_versions
    }
}
