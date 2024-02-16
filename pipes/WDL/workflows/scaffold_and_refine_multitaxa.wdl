version 1.0

import "../tasks/tasks_assembly.wdl" as assembly
import "../tasks/tasks_metagenomics.wdl" as metagenomics
import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_utils.wdl" as utils
import "assemble_refbased.wdl" as assemble_refbased

workflow scaffold_and_refine_multitaxa {
    meta {
        description: "Scaffold de novo contigs against a set of possible references and subsequently polish with reads."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        String  sample_id
        File    reads_unmapped_bam

        File    taxid_to_ref_accessions_tsv
        File?   focal_report_tsv
        File?   ncbi_taxdump_tgz

        # Float    min_pct_reference_covered = 0.1
    }

    Int  min_scaffold_unambig = 10

    # if kraken reports are available, filter scaffold list to observed hits (output might be empty!)
    if(defined(focal_report_tsv) && defined(ncbi_taxdump_tgz)) {
        call metagenomics.filter_refs_to_found_taxa {
            input:
                taxid_to_ref_accessions_tsv = taxid_to_ref_accessions_tsv,
                taxdump_tgz = select_first([ncbi_taxdump_tgz]),
                focal_report_tsv = select_first([focal_report_tsv])
        }
    }

    Array[Array[String]] taxid_to_ref_accessions = read_tsv(select_first([filter_refs_to_found_taxa.filtered_taxid_to_ref_accessions_tsv, taxid_to_ref_accessions_tsv]))
    Array[String] assembly_header = ["sample_id", "taxid", "tax_name", "assembly_fasta", "aligned_only_reads_bam", "coverage_plot", "assembly_length", "assembly_length_unambiguous", "reads_aligned", "mean_coverage", "percent_reference_covered", "intermediate_gapfill_fasta", "assembly_preimpute_length_unambiguous", "replicate_concordant_sites", "replicate_discordant_snps", "replicate_discordant_indels", "replicate_discordant_vcf", "isnvsFile", "aligned_bam", "coverage_tsv", "read_pairs_aligned", "bases_aligned"]

    scatter(taxon in taxid_to_ref_accessions) {
        # taxon = [taxid, taxname, semicolon_delim_accession_list]
        if(length(taxon)>1) { # <-- workaround for serious bug in cromwell's read_tsv on empty files
            call utils.string_split {
                input:
                    joined_string = taxon[2],
                    delimiter = ";"
            }
            call ncbi.download_annotations {
                input:
                    accessions = string_split.tokens,
                    combined_out_prefix = taxon[0]
            }
            call assembly.scaffold {
                input:
                    reads_bam = reads_unmapped_bam,
                    reference_genome_fasta = [download_annotations.combined_fasta],
                    min_length_fraction = 0,
                    min_unambig = 0,
                    allow_incomplete_output = true
            }
            if (scaffold.assembly_preimpute_length_unambiguous > min_scaffold_unambig) {
                # polish de novo assembly with reads
                call assemble_refbased.assemble_refbased as refine {
                    input:
                        reads_unmapped_bams = [reads_unmapped_bam],
                        reference_fasta     = scaffold.scaffold_fasta,
                        sample_name         = sample_id
                }
                String assembly_method_denovo = "viral-ngs/assemble_denovo"
            }
            if (scaffold.assembly_preimpute_length_unambiguous <= min_scaffold_unambig) {
                # fall back to refbased assembly if de novo fails
                call assemble_refbased.assemble_refbased as ref_based {
                    input:
                        reads_unmapped_bams = [reads_unmapped_bam],
                        reference_fasta     = download_annotations.combined_fasta,
                        sample_name         = sample_id
                }
                String assembly_method_refbased = "viral-ngs/assemble_refbased"
            }

            # TO DO: if percent_reference_covered > some threshold, run ncbi.rename_fasta_header and ncbi.align_and_annot_transfer_single
            # TO DO: if biosample attributes file provided, run ncbi.biosample_to_genbank

            String taxid = taxon[0]
            String tax_name = taxon[1]
            Int    assembly_length_unambiguous = select_first([refine.assembly_length_unambiguous, ref_based.assembly_length_unambiguous])
            Float  percent_reference_covered = 1.0 * assembly_length_unambiguous / scaffold.reference_length
            File   assembly_fasta = select_first([refine.assembly_fasta, ref_based.assembly_fasta])
            Map[String, String] stats_by_taxon = {
                "sample_id" : sample_id,
                "taxid" : taxid,
                "tax_name" : tax_name,

                "assembly_fasta" :              assembly_fasta,
                "aligned_only_reads_bam" :      select_first([refine.align_to_self_merged_aligned_only_bam, ref_based.align_to_self_merged_aligned_only_bam]),
                "coverage_plot" :               select_first([refine.align_to_self_merged_coverage_plot, ref_based.align_to_self_merged_coverage_plot]),
                "assembly_length" :             select_first([refine.assembly_length, ref_based.assembly_length]),
                "assembly_length_unambiguous" : assembly_length_unambiguous,
                "reads_aligned" :               select_first([refine.align_to_self_merged_reads_aligned, ref_based.align_to_self_merged_reads_aligned]),
                "mean_coverage" :               select_first([refine.align_to_self_merged_mean_coverage, ref_based.align_to_self_merged_mean_coverage]),
                "percent_reference_covered" :   percent_reference_covered,

                "intermediate_gapfill_fasta" :            scaffold.intermediate_gapfill_fasta,
                "assembly_preimpute_length_unambiguous" : scaffold.assembly_preimpute_length_unambiguous,

                "replicate_concordant_sites" :  select_first([refine.replicate_concordant_sites, ref_based.replicate_concordant_sites]),
                "replicate_discordant_snps" :   select_first([refine.replicate_discordant_snps, ref_based.replicate_discordant_snps]),
                "replicate_discordant_indels" : select_first([refine.replicate_discordant_indels, ref_based.replicate_discordant_indels]),
                "replicate_discordant_vcf" :    select_first([refine.replicate_discordant_vcf, ref_based.replicate_discordant_vcf]),

                "isnvsFile" :          select_first([refine.align_to_self_isnvs_vcf, ref_based.align_to_self_isnvs_vcf]),
                "aligned_bam" :        select_first([refine.align_to_self_merged_aligned_only_bam, ref_based.align_to_self_merged_aligned_only_bam]),
                "coverage_tsv" :       select_first([refine.align_to_self_merged_coverage_tsv, ref_based.align_to_self_merged_coverage_tsv]),
                "read_pairs_aligned" : select_first([refine.align_to_self_merged_read_pairs_aligned, ref_based.align_to_self_merged_read_pairs_aligned]),
                "bases_aligned" :      select_first([refine.align_to_self_merged_bases_aligned, ref_based.align_to_self_merged_bases_aligned]),

                "assembly_method" :    select_first([assembly_method_denovo, assembly_method_refbased])
            }

            scatter(h in assembly_header) {
                String stat_by_taxon = stats_by_taxon[h]
            }
        }
    }
 
    ### summary stats
    call utils.concatenate {
      input:
        infiles     = [write_tsv([assembly_header]), write_tsv(select_all(stat_by_taxon))],
        output_name = "assembly_metadata-~{sample_id}.tsv"
    }

    output {
        Array[Map[String,String]] assembly_stats_by_taxon  = select_all(stats_by_taxon)
        File   assembly_stats_by_taxon_tsv                 = concatenate.combined
        String assembly_method                             = "viral-ngs/scaffold_and_refine_multitaxa"

        ## TO DO: some summary stats on stats_by_taxon: how many rows, numbers from the best row, etc
        #String assembly_top_taxon_id               = 
        #String assembly_top_length_unambiguous     = 
        #Float  assembly_top_pct_ref_cov            = 
        #File   assembly_top_fasta                  = 
        Array[String] assembly_all_taxids          = select_all(taxid)
        Array[String] assembly_all_taxnames        = select_all(tax_name)
        Array[Int]    assembly_all_lengths_unambig = select_all(assembly_length_unambiguous)
        Array[Float]  assembly_all_pct_ref_cov     = select_all(percent_reference_covered)
        Array[File]   assembly_all_fastas          = select_all(assembly_fasta)
    }
}
