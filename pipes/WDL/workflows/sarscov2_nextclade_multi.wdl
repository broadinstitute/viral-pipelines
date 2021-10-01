version 1.0

import "../tasks/tasks_sarscov2.wdl" as sarscov2

workflow sarscov2_nextclade_multi {
    meta {
        description: "Create Nextclade visualizations on many SARS-CoV-2 genomes"
    }

    input {
        Array[File]+ genome_fastas
    }

    call sarscov2.nextclade_many_samples {
        input:
            genome_fastas = genome_fastas,
            dataset_name  = "sars-cov-2"
    }

    output {
        File nextclade_tsv  = nextclade_many_samples.nextclade_tsv
        File nextclade_json = nextclade_many_samples.nextclade_json
        File auspice_json   = nextclade_many_samples.auspice_json
    }
}
