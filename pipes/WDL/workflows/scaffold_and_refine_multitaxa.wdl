version 1.0

import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_utils.wdl" as utils
import "assemble_refbased.wdl" as assemble_refbased

workflow scaffold_and_refine_multitaxa {
    meta {
        description: "Scaffold de novo contigs against a set of possible references and subsequently polish with reads. This workflow accepts a very large set of input reference genomes. It will subset the reference genomes to those with ANI hits to the provided contigs/MAGs and cluster the reference hits by any ANI similarity to each other. It will choose the top reference from each cluster and produce one assembly for each cluster. This is intended to allow for the presence of multiple diverse viral taxa (coinfections) while forcing a choice of the best assembly from groups of related reference genomes."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        String  sample_id
        String? sample_name
        File    reads_unmapped_bam
        File    contigs_fasta

        File    taxid_to_ref_accessions_tsv
    }

    Int    min_scaffold_unambig = 10 # in base-pairs; any scaffolded assembly < this length will not be refined/polished
    String sample_original_name = select_first([sample_name, sample_id])

    # download (multi-segment) genomes for each reference, fasta filename = colon-concatenated accession list
    scatter(taxon in read_tsv(taxid_to_ref_accessions_tsv)) {
        # taxon = [taxid, taxname, semicolon_delim_accession_list]
        call utils.string_split {
            input:
                joined_string = taxon[2],
                delimiter = ":"
        }
        call ncbi.download_annotations {
            input:
                accessions = string_split.tokens,
                combined_out_prefix = taxon[2]
        }
    }

    # subset reference genomes to those with ANI hits to contigs and cluster reference hits by any ANI similarity to each other
    call assembly.select_references {
        input:
            reference_genomes_fastas = download_annotations.combined_fasta,
            contigs_fasta = contigs_fasta
    }

    # assemble and produce stats for every reference cluster
    Array[String] assembly_header = ["entity:assembly_id", "assembly_name", "sample_id", "sample_name", "taxid", "tax_name", "assembly_fasta", "aligned_only_reads_bam", "coverage_plot", "assembly_length", "assembly_length_unambiguous", "reads_aligned", "mean_coverage", "percent_reference_covered", "scaffolding_num_segments_recovered", "reference_num_segments_required", "reference_length", "reference_accessions", "skani_num_ref_clusters", "skani_this_cluster_num_refs", "skani_dist_tsv", "scaffolding_ani", "scaffolding_pct_ref_cov", "intermediate_gapfill_fasta", "assembly_preimpute_length_unambiguous", "replicate_concordant_sites", "replicate_discordant_snps", "replicate_discordant_indels", "replicate_discordant_vcf", "isnvsFile", "aligned_bam", "coverage_tsv", "read_pairs_aligned", "bases_aligned", "coverage_genbank", "assembly_method", "sample"]
    scatter(ref_cluster_tar in select_references.matched_reference_clusters_fastas_tars) {

        call utils.tar_extract {
            input:
                tar_file = ref_cluster_tar
        }

        # assemble (scaffold-and-refine) genome against this reference cluster
        call assembly.scaffold {
            input:
                reads_bam = reads_unmapped_bam,
                contigs_fasta = contigs_fasta,
                reference_genome_fasta = tar_extract.files,
                min_length_fraction = 0,
                min_unambig = 0,
                allow_incomplete_output = true
        }
        if (scaffold.assembly_preimpute_length_unambiguous > min_scaffold_unambig) {
            # polish de novo assembly with reads
            call assemble_refbased.assemble_refbased as refine {
                input:
                    reads_unmapped_bams  = [reads_unmapped_bam],
                    reference_fasta      = scaffold.scaffold_fasta,
                    sample_name          = sample_id,
                    sample_original_name = sample_original_name
            }
            call reports.coverage_report as coverage_self {
                input:
                    mapped_bams = [refine.align_to_self_merged_aligned_only_bam],
                    mapped_bam_idx = []
            }
            call utils.tsv_drop_cols as coverage_two_col {
                input:
                    in_tsv = coverage_self.coverage_report,
                    drop_cols = ["aln2self_cov_median", "aln2self_cov_mean_non0", "aln2self_cov_1X", "aln2self_cov_5X", "aln2self_cov_20X", "aln2self_cov_100X"]
            }
        }

        # get taxid and taxname from taxid_to_ref_accessions_tsv
        call utils.fetch_row_from_tsv as tax_lookup {
            input:
                tsv = taxid_to_ref_accessions_tsv,
                idx_col = "accessions",
                idx_val = scaffold.scaffolding_chosen_ref_basename,
                add_header = ["taxid", "taxname", "accessions"]
        }
        String taxid = tax_lookup.map["taxid"]
        String tax_name = tax_lookup.map["taxname"]

        # build output tsv row
        Int    assembly_length_unambiguous = select_first([refine.assembly_length_unambiguous, 0])
        Float  percent_reference_covered = 1.0 * assembly_length_unambiguous / scaffold.reference_length
        File   assembly_fasta = select_first([refine.assembly_fasta, scaffold.intermediate_gapfill_fasta])
        Map[String, String] stats_by_taxon = {
            "entity:assembly_id" : sample_id + "-" + taxid,
            "assembly_name" :      tax_name + ": " + sample_original_name,
            "sample_id" :          sample_id,
            "sample_name" :        sample_original_name,
            "taxid" :              taxid,
            "tax_name" :           tax_name,

            "assembly_fasta" :              assembly_fasta,
            "aligned_only_reads_bam" :      select_first([refine.align_to_self_merged_aligned_only_bam, ""]),
            "coverage_plot" :               select_first([refine.align_to_self_merged_coverage_plot, ""]),
            "assembly_length" :             select_first([refine.assembly_length, "0"]),
            "assembly_length_unambiguous" : assembly_length_unambiguous,
            "reads_aligned" :               select_first([refine.align_to_self_merged_reads_aligned, "0"]),
            "mean_coverage" :               select_first([refine.align_to_self_merged_mean_coverage, "0"]),
            "percent_reference_covered" :   percent_reference_covered,
            "scaffolding_num_segments_recovered" : scaffold.assembly_num_segments_recovered,
            "reference_num_segments_required" : scaffold.reference_num_segments_required,
            "reference_length" :            scaffold.reference_length,
            "reference_accessions" :        tax_lookup.map["accessions"],

            "skani_num_ref_clusters" :      length(select_references.matched_reference_clusters_fastas_tars),
            "skani_this_cluster_num_refs" : length(tar_extract.files),
            "skani_dist_tsv" :              scaffold.scaffolding_stats,
            "scaffolding_ani" :             scaffold.skani_ani,
            "scaffolding_pct_ref_cov" :     scaffold.skani_ref_aligned_frac,

            "intermediate_gapfill_fasta" :            scaffold.intermediate_gapfill_fasta,
            "assembly_preimpute_length_unambiguous" : scaffold.assembly_preimpute_length_unambiguous,

            "replicate_concordant_sites" :  select_first([refine.replicate_concordant_sites, "0"]),
            "replicate_discordant_snps" :   select_first([refine.replicate_discordant_snps, "0"]),
            "replicate_discordant_indels" : select_first([refine.replicate_discordant_indels, "0"]),
            "replicate_discordant_vcf" :    select_first([refine.replicate_discordant_vcf, ""]),

            "isnvsFile" :          select_first([refine.align_to_self_isnvs_vcf, ""]),
            "aligned_bam" :        select_first([refine.align_to_self_merged_aligned_only_bam, ""]),
            "coverage_tsv" :       select_first([refine.align_to_self_merged_coverage_tsv, ""]),
            "read_pairs_aligned" : select_first([refine.align_to_self_merged_read_pairs_aligned, "0"]),
            "bases_aligned" :      select_first([refine.align_to_self_merged_bases_aligned, "0"]),
            "coverage_genbank" :   select_first([coverage_two_col.out_tsv, ""]),

            "assembly_method" :    "viral-ngs/assemble_denovo",

            "sample":              '{"entityType":"sample","entityName":"' + sample_id + '"}'
        }

        if(assembly_length_unambiguous > min_scaffold_unambig) {
            scatter(h in assembly_header) {
                String stat_by_taxon = stats_by_taxon[h]
            }
        }
    }

    ### summary stats
    if (length(select_all(stat_by_taxon)) > 0) {
        call utils.concatenate as assembly_stats_non_empty {
            input:
                infiles     = [write_tsv([assembly_header]), write_tsv(select_all(stat_by_taxon))],
                output_name = "assembly_metadata-~{sample_id}.tsv"
        }
    }
    if (length(select_all(stat_by_taxon)) == 0) {
        call utils.concatenate as assembly_stats_empty {
            input:
                infiles     = [write_tsv([assembly_header])],
                output_name = "assembly_metadata-~{sample_id}.tsv"
        }
    }

    output {
        Array[Map[String,String]] assembly_stats_by_taxon  = stats_by_taxon
        File   assembly_stats_by_taxon_tsv                 = select_first([assembly_stats_non_empty.combined, assembly_stats_empty.combined])
        String assembly_method                             = "viral-ngs/scaffold_and_refine_multitaxa"

        #String assembly_top_taxon_id               = select_references.top_matches_per_cluster_basenames[0]
        Int    skani_num_ref_clusters              = length(select_references.matched_reference_clusters_fastas_tars)
        File   skani_contigs_to_refs_dist_tsv      = select_references.skani_dist_full_tsv

        Array[String] assembly_all_taxids          = taxid
        Array[String] assembly_all_taxnames        = tax_name
        Array[Int]    assembly_all_lengths_unambig = assembly_length_unambiguous
        Array[Float]  assembly_all_pct_ref_cov     = percent_reference_covered
        Array[File]   assembly_all_fastas          = assembly_fasta
    }
}
