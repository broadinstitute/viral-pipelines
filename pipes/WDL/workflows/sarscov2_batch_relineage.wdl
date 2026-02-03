version 1.1

import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_sarscov2.wdl" as sarscov2
import "../tasks/tasks_utils.wdl" as utils

workflow sarscov2_batch_relineage {
    meta {
        description: "Re-call Nextclade and Pangolin lineages on a flowcell's worth of SARS-CoV-2 genomes"
        allowNestedInputs: true
    }

    input {
        String      flowcell_id
        Array[File] genomes_fasta
        File        metadata_annotated_tsv
        File        metadata_raw_tsv
        Int         min_genome_bases = 24000
    }

    call utils.concatenate {
        input:
            infiles     = genomes_fasta,
            output_name = "all-genomes.fasta",
            cpus        = 16
    }

    call utils.filter_sequences_by_length {
        input:
            sequences_fasta = concatenate.combined,
            min_non_N       = min_genome_bases
    }

    call utils.fasta_to_ids {
        input:
            sequences_fasta = filter_sequences_by_length.filtered_fasta
    }

    call nextstrain.nextclade_many_samples {
        input:
            genome_fastas               = [filter_sequences_by_length.filtered_fasta],
            genome_ids_setdefault_blank = fasta_to_ids.ids_txt,
            basename                    = "nextclade-~{flowcell_id}",
            dataset_name                = "sars-cov-2"
    }

    call sarscov2.pangolin_many_samples {
        input:
            genome_fastas = [filter_sequences_by_length.filtered_fasta],
            basename      = "pangolin-~{flowcell_id}"
    }

    scatter(sample_sanitized in read_lines(fasta_to_ids.ids_txt)) {
        Array[String] metadata = [
            sample_sanitized,
            pangolin_many_samples.pango_lineage[sample_sanitized],
            pangolin_many_samples.scorpio_call[sample_sanitized],
            nextclade_many_samples.nextclade_clade[sample_sanitized],
            nextclade_many_samples.aa_subs_csv[sample_sanitized],
            nextclade_many_samples.aa_dels_csv[sample_sanitized],
            pangolin_many_samples.pangolin_versions,
            nextclade_many_samples.nextclade_version
        ]
    }
    Array[String] meta_header = [
        'sample_sanitized',
        'pango_lineage', 'scorpio_call',
        'nextclade_clade', 'nextclade_aa_subs', 'nextclade_aa_dels',
        'pangolin_version', 'nextclade_version'
    ]

    call utils.today

    call utils.tsv_join as merge_raw {
        input:
            input_tsvs = [write_tsv(flatten([[meta_header], metadata])),
                metadata_raw_tsv],
            id_col = "sample_sanitized",
            out_suffix = ".tsv",
            out_basename = basename(metadata_raw_tsv, '.tsv') + ".relineage_~{today.date}"
    }

    call utils.tsv_join as merge_annotated {
        input:
            input_tsvs = [write_tsv(flatten([[meta_header], metadata])),
                metadata_annotated_tsv],
            id_col = "sample_sanitized",
            out_suffix = ".tsv",
            out_basename = basename(metadata_annotated_tsv, '.tsv') + ".relineage_~{today.date}"
    }

    output {
        File   assembly_stats_relineage_tsv = merge_raw.out_tsv
        File   assembly_stats_final_relineage_tsv = merge_annotated.out_tsv
        File   nextclade_all_json           = nextclade_many_samples.nextclade_json
        File   nextclade_all_tsv            = nextclade_many_samples.nextclade_tsv
        File   nextclade_auspice_json       = nextclade_many_samples.auspice_json
        File   nextalign_msa                = nextclade_many_samples.nextalign_msa
        File   pangolin_report              = pangolin_many_samples.pango_lineage_report
        File   pangolin_msa                 = pangolin_many_samples.msa_fasta
    }
}
