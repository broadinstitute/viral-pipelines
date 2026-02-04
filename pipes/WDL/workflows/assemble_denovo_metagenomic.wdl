version 1.0

import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_read_utils.wdl" as read_utils
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_utils.wdl" as utils
import "assemble_refbased.wdl" as assemble_refbased

workflow assemble_denovo_metagenomic {
    meta {
        description: "Performs viral de novo assembly on metagenomic reads against a large range of possible reference genomes. Runs raw reads through taxonomic classification (Kraken2), human read depletion (based on Kraken2), de novo assembly (SPAdes), and FASTQC/multiQC of reads. Scaffold de novo contigs against a set of possible references and subsequently polish with reads. This workflow accepts a very large set of input reference genomes. It will subset the reference genomes to those with ANI hits to the provided contigs/MAGs and cluster the reference hits by any ANI similarity to each other. It will choose the top reference from each cluster and produce one assembly for each cluster. This is intended to allow for the presence of multiple diverse viral taxa (coinfections) while forcing a choice of the best assembly from groups of related reference genomes."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        String        sample_id
        String?       sample_name
        String?       biosample_accession

        Array[File]+  reads_bams

        Array[String] batch_id_list

        File          ncbi_taxdump_tgz

        File?         spikein_db
        File          trim_clip_db

        File          kraken2_db_tgz
        File          krona_taxonomy_db_kraken2_tgz

        File          taxid_to_ref_accessions_tsv

        Array[String] taxa_to_dehost         = ["Vertebrata"]
        Array[String] taxa_to_avoid_assembly = ["Vertebrata", "other sequences", "Bacteria"]

        Int           min_reads_for_rmdup    =  5000000
        Int           max_reads_for_assembly = 10000000

        String        table_name = "sample"
    }

    Int    min_scaffold_unambig = 300 # in base-pairs; any scaffolded assembly < this length will not be refined/polished
    String sample_original_name = select_first([sample_name, sample_id])

    parameter_meta {
        reads_bams: {
          description: "Reads to classify. May be unmapped or mapped or both, paired-end or single-end. Multiple input files will be merged first.",
          patterns: ["*.bam"]
        }
        spikein_db: {
          description: "ERCC spike-in sequences",
          patterns: ["*.fasta", "*.fasta.gz", "*.fasta.zst"]
        }
        trim_clip_db: {
          description: "Adapter sequences to remove via trimmomatic prior to SPAdes assembly",
          patterns: ["*.fasta", "*.fasta.gz", "*.fasta.zst"]
        }
        kraken2_db_tgz: {
          description: "Pre-built Kraken database tarball containing three files: hash.k2d, opts.k2d, and taxo.k2d.",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        krona_taxonomy_db_kraken2_tgz: {
          description: "Krona taxonomy database containing a single file: taxonomy.tab, or possibly just a compressed taxonomy.tab",
          patterns: ["*.tab.zst", "*.tab.gz", "*.tab", "*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        ncbi_taxdump_tgz: {
          description: "An NCBI taxdump.tar.gz file that contains, at the minimum, a nodes.dmp and names.dmp file.",
          patterns: ["*.tar.gz", "*.tar.lz4", "*.tar.bz2", "*.tar.zst"]
        }
        cleaned_fastqc: { 
          description: "Output cleaned fastqc reports in HTML.",
          category: "other"
        }
        kraken2_krona_plot: {
          description:"Visualize the results of the Kraken2 analysis with Krona, which disaplys taxonmic hierarchiral data in multi-layerd pie.",
          category:"other"
        }
        kraken2_summary_report: {
          description: "Kraken report output file.",
          category: "other"
        }
        raw_fastqc:{
          description: "Merged raw fastqc reads.",
          category: "other"
        }

    }

    call utils.unique_strings as unique_batch_ids {
        input:
            strings = batch_id_list
    }

    if(length(reads_bams) > 1) {
      call read_utils.merge_and_reheader_bams as merge_raw_reads {
        input:
            in_bams      = reads_bams,
            out_basename = sample_id
      }
    }
    File reads_bam = select_first([merge_raw_reads.out_bam, reads_bams[0]])

    if(defined(spikein_db)) {
        call reports.align_and_count as spikein {
            input:
                reads_bam = reads_bam,
                ref_db    = select_first([spikein_db])
        }
    }
    call metagenomics.kraken2 as kraken2 {
        input:
            reads_bam             = reads_bam,
            kraken2_db_tgz        = kraken2_db_tgz,
            krona_taxonomy_db_tgz = krona_taxonomy_db_kraken2_tgz
    }
    call metagenomics.filter_bam_to_taxa as deplete {
        input:
            classified_bam          = reads_bam,
            classified_reads_txt_gz = kraken2.kraken2_reads_report,
            ncbi_taxonomy_db_tgz    = ncbi_taxdump_tgz,
            exclude_taxa            = true,
            taxonomic_names         = taxa_to_dehost,
            out_filename_suffix     = "hs_depleted"
    }
    call metagenomics.filter_bam_to_taxa as filter_acellular {
        input:
            classified_bam          = reads_bam,
            classified_reads_txt_gz = kraken2.kraken2_reads_report,
            ncbi_taxonomy_db_tgz    = ncbi_taxdump_tgz,
            exclude_taxa            = true,
            taxonomic_names         = taxa_to_avoid_assembly,
            out_filename_suffix     = "acellular"
    }
    call read_utils.bbnorm_bam {
        input:
            reads_bam        = filter_acellular.bam_filtered_to_taxa,
            min_input_reads  = min_reads_for_rmdup,
            max_output_reads = max_reads_for_assembly
    }

    call assembly.assemble as spades {
        input:
            reads_unmapped_bam = bbnorm_bam.bbnorm_bam,
            trim_clip_db       = trim_clip_db,
            always_succeed     = true
    }
    call metagenomics.report_primary_kraken_taxa {
        input:
            kraken_summary_report = kraken2.kraken2_summary_report
    }

    # download (multi-segment) genomes for each reference, fasta filename = dash-concatenated accession list
    call ncbi.download_ref_genomes_from_tsv {
        input:
            ref_genomes_tsv = taxid_to_ref_accessions_tsv
    }

    # subset reference genomes to those with ANI hits to contigs and cluster reference hits by any ANI similarity to each other
    call assembly.select_references {
        input:
            reference_genomes_fastas = download_ref_genomes_from_tsv.ref_genomes_fastas,
            contigs_fasta = spades.contigs_fasta
    }

    # assemble and produce stats for every reference cluster
    Array[String] assembly_header = ["entity:assembly_id", "assembly_name", "sample_id", "sample_name", "taxid", "tax_name", "tax_shortname", "assembly_fasta", "aligned_only_reads_bam", "coverage_plot", "assembly_length", "assembly_length_unambiguous", "reads_aligned", "mean_coverage", "percent_reference_covered", "scaffolding_num_segments_recovered", "reference_num_segments_required", "reference_length", "reference_accessions", "skani_num_ref_clusters", "skani_this_cluster_num_refs", "skani_dist_tsv", "scaffolding_ani", "scaffolding_pct_ref_cov", "intermediate_gapfill_fasta", "assembly_preimpute_length_unambiguous", "replicate_concordant_sites", "replicate_discordant_snps", "replicate_discordant_indels", "replicate_discordant_vcf", "isnvsFile", "aligned_bam", "coverage_tsv", "read_pairs_aligned", "bases_aligned", "assembly_method", "assembly_method_version", "biosample_accession", "batch_ids", "~{table_name}"]
    scatter(ref_cluster_tar in select_references.matched_reference_clusters_fastas_tars) {

        call utils.tar_extract {
            input:
                tar_file = ref_cluster_tar
        }

        # assemble (scaffold-and-refine) genome against this reference cluster
        call assembly.scaffold {
            input:
                reads_bam = bbnorm_bam.bbnorm_bam,
                contigs_fasta = spades.contigs_fasta,
                reference_genome_fasta = tar_extract.files,
                min_length_fraction = 0,
                min_unambig = 0,
                allow_incomplete_output = true
        }

        # get taxid and taxname from taxid_to_ref_accessions_tsv
        call utils.fetch_row_from_tsv as tax_lookup {
            input:
                tsv = taxid_to_ref_accessions_tsv,
                idx_col = "accessions",
                idx_val = sub(scaffold.scaffolding_chosen_ref_basename, "-", ":"),
                add_header = ["tax_id", "isolate_prefix", "taxname", "accessions"]
        }
        String taxid = tax_lookup.map["tax_id"]
        String tax_name = tax_lookup.map["taxname"]
        String isolate_prefix = tax_lookup.map["isolate_prefix"]

        # polish de novo assembly with reads
        if (scaffold.assembly_preimpute_length_unambiguous > min_scaffold_unambig) {
            call assemble_refbased.assemble_refbased as refine {
                input:
                    reads_unmapped_bams  = [filter_acellular.bam_filtered_to_taxa],
                    reference_fasta      = scaffold.scaffold_fasta,
                    sample_name          = sample_id + "-" + taxid
            }
        }

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
            "tax_shortname" :      isolate_prefix,

            "assembly_fasta" :              "~{assembly_fasta}",
            "aligned_only_reads_bam" :      "~{select_first([refine.align_to_self_merged_aligned_only_bam, ''])}",
            "coverage_plot" :               "~{select_first([refine.align_to_self_merged_coverage_plot, ''])}",
            "assembly_length" :             "~{select_first([refine.assembly_length, 0])}",
            "assembly_length_unambiguous" : "~{assembly_length_unambiguous}",
            "reads_aligned" :               "~{select_first([refine.align_to_self_merged_reads_aligned, 0])}",
            "mean_coverage" :               "~{select_first([refine.align_to_self_merged_mean_coverage, 0])}",
            "percent_reference_covered" :   "~{percent_reference_covered}",
            "scaffolding_num_segments_recovered" : "~{scaffold.assembly_num_segments_recovered}",
            "reference_num_segments_required" : "~{scaffold.reference_num_segments_required}",
            "reference_length" :            "~{scaffold.reference_length}",
            "reference_accessions" :        tax_lookup.map["accessions"],

            "skani_num_ref_clusters" :      "~{length(select_references.matched_reference_clusters_fastas_tars)}",
            "skani_this_cluster_num_refs" : "~{length(tar_extract.files)}",
            "skani_dist_tsv" :              "~{scaffold.scaffolding_stats}",
            "scaffolding_ani" :             "~{scaffold.skani_ani}",
            "scaffolding_pct_ref_cov" :     "~{scaffold.skani_ref_aligned_frac}",

            "intermediate_gapfill_fasta" :            "~{scaffold.intermediate_gapfill_fasta}",
            "assembly_preimpute_length_unambiguous" : "~{scaffold.assembly_preimpute_length_unambiguous}",

            "replicate_concordant_sites" :  "~{select_first([refine.replicate_concordant_sites, 0])}",
            "replicate_discordant_snps" :   "~{select_first([refine.replicate_discordant_snps, 0])}",
            "replicate_discordant_indels" : "~{select_first([refine.replicate_discordant_indels, 0])}",
            "replicate_discordant_vcf" :    "~{select_first([refine.replicate_discordant_vcf, ''])}",

            "isnvsFile" :          "~{select_first([refine.align_to_self_isnvs_vcf, ''])}",
            "aligned_bam" :        "~{select_first([refine.align_to_self_merged_aligned_only_bam, ''])}",
            "coverage_tsv" :       "~{select_first([refine.align_to_self_merged_coverage_tsv, ''])}",
            "read_pairs_aligned" : "~{select_first([refine.align_to_self_merged_read_pairs_aligned, 0])}",
            "bases_aligned" :      "~{select_first([refine.align_to_self_merged_bases_aligned, 0])}",

            "assembly_method" :         "viral-ngs/assemble_denovo",
            "assembly_method_version" : scaffold.viralngs_version,

            "biosample_accession" :     "~{select_first([biosample_accession, ''])}",

            "batch_ids" :          unique_batch_ids.sorted_unique_joined,

            "~{table_name}":            '{"entityType":"~{table_name}","entityName":"' + sample_id + '"}'
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
        File   reads_dehosted_ubam             = deplete.bam_filtered_to_taxa
        File   reads_acellular_ubam            = filter_acellular.bam_filtered_to_taxa
        File   reads_assembly_input_ubam       = bbnorm_bam.bbnorm_bam
        File   contigs_fasta                   = spades.contigs_fasta

        Int    read_counts_raw                 = deplete.classified_taxonomic_filter_read_count_pre
        Int    read_counts_dehosted            = deplete.classified_taxonomic_filter_read_count_post
        Int    read_counts_acellular           = filter_acellular.classified_taxonomic_filter_read_count_post
        Int    read_counts_assembly_input      = bbnorm_bam.bbnorm_read_count_post
        Int    read_counts_prespades_subsample = spades.subsample_read_count
        
        File   kraken2_summary_report          = kraken2.kraken2_summary_report
        File   kraken2_reads_report            = kraken2.kraken2_reads_report
        File   kraken2_krona_plot              = kraken2.krona_report_html
        File   kraken2_top_taxa_report         = report_primary_kraken_taxa.ranked_focal_report
        String kraken2_focal_taxon_name        = report_primary_kraken_taxa.focal_tax_name
        Int    kraken2_focal_total_reads       = report_primary_kraken_taxa.total_focal_reads
        String kraken2_top_taxon_id            = report_primary_kraken_taxa.tax_id
        String kraken2_top_taxon_name          = report_primary_kraken_taxa.tax_name
        String kraken2_top_taxon_rank          = report_primary_kraken_taxa.tax_rank
        Int    kraken2_top_taxon_num_reads     = report_primary_kraken_taxa.num_reads
        Float  kraken2_top_taxon_pct_of_focal  = report_primary_kraken_taxa.percent_of_focal

        File?   spikein_report                 = spikein.report
        String? spikein_tophit                 = spikein.top_hit_id
        String? spikein_pct_of_total_reads     = spikein.pct_total_reads_mapped
        String? spikein_pct_lesser_hits        = spikein.pct_lesser_hits_of_mapped
        
        String viral_classify_version          = kraken2.viralngs_version
        String viral_assemble_version          = spades.viralngs_version

        Array[Map[String,String]] assembly_stats_by_taxon  = stats_by_taxon
        File   assembly_stats_by_taxon_tsv                 = select_first([assembly_stats_non_empty.combined, assembly_stats_empty.combined])
        String assembly_method                             = "viral-ngs/scaffold_and_refine_multitaxa"

        Int    skani_num_ref_clusters              = length(select_references.matched_reference_clusters_fastas_tars)
        File   skani_contigs_to_refs_dist_tsv      = select_references.skani_dist_full_tsv

        Array[String] assembly_all_taxids          = taxid
        Array[String] assembly_all_taxnames        = tax_name
        Array[Int]    assembly_all_lengths_unambig = assembly_length_unambiguous
        Array[Float]  assembly_all_pct_ref_cov     = percent_reference_covered
        Array[File]   assembly_all_fastas          = assembly_fasta

        String  batch_ids                          = unique_batch_ids.sorted_unique_joined
    }
}
