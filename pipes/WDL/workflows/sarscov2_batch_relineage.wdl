version 1.0

import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_sarscov2.wdl" as sarscov2
import "../tasks/tasks_utils.wdl" as utils

workflow sarscov2_batch_relineage {
    meta {
        description: "Re-call Nextclade and Pangolin lineages on a flowcell's worth of SARS-CoV-2 genomes"
    }

    input {
        Array[File] genomes_fasta
        File        metadata_annotated_tsv
        File        metadata_raw_tsv
        Int         min_genome_bases = 24000
    }

    scatter(genome_fasta in genomes_fasta) {
        call reports.assembly_bases {
            input: fasta = genome_fasta
        }

        if (assembly_bases.assembly_length_unambiguous >= min_genome_bases) {
            call sarscov2.nextclade_one_sample {
                input:
                    genome_fasta = genome_fasta
            }
            call sarscov2.pangolin_one_sample {
                input:
                    genome_fasta = genome_fasta
            }
            Array[String] metadata = [
                basename(genome_fasta, '.fasta'),
                pangolin_one_sample.pango_lineage,
                nextclade_one_sample.nextclade_clade,
                nextclade_one_sample.aa_subs_csv,
                nextclade_one_sample.aa_dels_csv,
                pangolin_one_sample.version,
                nextclade_one_sample.nextclade_tsv,#
                nextclade_one_sample.nextclade_json,#
                pangolin_one_sample.pango_lineage_report,#
            ]
        }
    }
    Array[String] meta_header = [
        'sample_sanitized',
        'pango_lineage', 'nextclade_clade', 'nextclade_aa_subs', 'nextclade_aa_dels', 'pangolin_version',
        'nextclade_tsv', 'nextclade_json', 'pangolin_csv'
    ]

    call utils.today

    call utils.tsv_join as merge_raw {
        input:
            input_tsvs = [write_tsv(flatten([[meta_header], select_all(metadata)])),
                metadata_raw_tsv],
            id_col = "sample_sanitized",
            out_basename = basename(metadata_raw_tsv, '.tsv') + "relineage_~{today.date}.tsv"
    }

    call utils.tsv_join as merge_annotated {
        input:
            input_tsvs = [write_tsv(flatten([[meta_header], select_all(metadata)])),
                metadata_annotated_tsv],
            id_col = "sample_sanitized",
            out_basename = basename(metadata_annotated_tsv, '.tsv') + "relineage_~{today.date}.tsv"
    }

    output {
        File   assembly_stats_relineage_tsv = merge_raw.out_tsv
        File   assembly_stats_final_relineage_tsv = merge_annotated.out_tsv
    }
}
