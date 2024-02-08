version 1.1

import "../tasks/tasks_assembly.wdl" as assembly
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

        #File    taxid_to_ref_accessions_json
        Map[Int,Array[String]+] taxid_to_ref_accessions = {
            208893:  ["KY654518.1"],  # RSV-A
            208895:  ["MZ516105.1"],  # RSV-B
            573824:  ["NC_038311.1"], # Rhino A1
            185900:  ["ON311191.1"],  # Rhino B27
            1418033: ["ON311169.1"],  # Rhino C15
            463676:  ["JN837686.2"],  # Rhino C45
            11137:   ["NC_002645.1"], # HCoV 229E
            290028:  ["NC_006577.2"], # HCoV HKU1
            277944:  ["NC_005831.2"], # HCoV NL63
            31631:   ["NC_006213.1"], # HCoV OC43
            2697049: ["NC_045512.2"], # SARS-CoV-2 Wuhan Hu-1
            641809:  ["NC_026438.1", "NC_026435.1", "NC_026437.1", "NC_026433.1", "NC_026436.1", "NC_026434.1", "NC_026431.1", "NC_026432.1"],  # Flu A/California/07/2009 H1N1
            335341:  ["NC_007373.1", "NC_007372.1", "NC_007371.1", "NC_007366.1", "NC_007369.1", "NC_007368.1", "NC_007367.1", "NC_007370.1"],  # Flu A/New York/392/2004 H3N2
            518987:  ["NC_002204.1", "NC_002205.1", "NC_002206.1", "NC_002207.1", "NC_002208.1", "NC_002209.1", "NC_002210.1", "NC_002211.1"],  # Flu B/Lee/1940
            162145:  ["NC_039199.1"], # metapneumo
            12730:   ["NC_003461.1"], # paraflu 1
            2560525: ["NC_003443.1"], # paraflu 2
            11216:   ["NC_001796.2"], # paraflu 3
            11224:   ["NC_021928.1"], # paraflu 4
            565302:  ["NC_011203.1"], # adenovirus B1
            565303:  ["NC_011202.1"], # adenovirus B2
            129951:  ["NC_001405.1"]  # adenovirus C
        }

        # Float    min_pct_reference_covered = 0.1
    }

    #Array[Pair[Int,Array[String]+]] taxid_to_ref_accessions = read_json(taxid_to_ref_accessions_json)
    Array[String] assembly_header = ["sample_id", "taxid", "assembly_fasta", "aligned_only_reads_bam", "coverage_plot", "assembly_length", "assembly_length_unambiguous", "reads_aligned", "mean_coverage", "percent_reference_covered", "intermediate_gapfill_fasta", "assembly_preimpute_length_unambiguous", "replicate_concordant_sites", "replicate_discordant_snps", "replicate_discordant_indels", "replicate_discordant_vcf", "isnvsFile", "aligned_bam", "coverage_tsv", "read_pairs_aligned", "bases_aligned"]

    scatter(taxid in keys(taxid_to_ref_accessions)) {
        call ncbi.download_annotations {
            input:
                accessions = taxid_to_ref_accessions[taxid],
                combined_out_prefix = taxid
        }
        call assembly.scaffold {
            input:
                reads_bam = reads_unmapped_bam,
                reference_genome_fasta = [download_annotations.combined_fasta],
                min_length_fraction = 0,
                min_unambig = 0,
                allow_incomplete_output = true
        }
        call assemble_refbased.assemble_refbased as refine {
            input:
                reads_unmapped_bams = [reads_unmapped_bam],
                reference_fasta     = scaffold.scaffold_fasta,
                sample_name         = sample_id
        }
        # to do: if percent_reference_covered > some threshold, run ncbi.rename_fasta_header and ncbi.align_and_annot_transfer_single
        # to do: if biosample attributes file provided, run ncbi.biosample_to_genbank

        if (refine.reference_genome_length > 0) {
            Float percent_reference_covered = 1.0 * refine.assembly_length_unambiguous / refine.reference_genome_length
        }

        Map[String, String] stats_by_taxon = {
            "sample_id" : sample_id,
            "taxid" : taxid,

            "assembly_fasta" : refine.assembly_fasta,
            "aligned_only_reads_bam" : refine.align_to_self_merged_aligned_only_bam,
            "coverage_plot" : refine.align_to_self_merged_coverage_plot,
            "assembly_length" : refine.assembly_length,
            "assembly_length_unambiguous" : refine.assembly_length_unambiguous,
            "reads_aligned" : refine.align_to_self_merged_reads_aligned,
            "mean_coverage" : refine.align_to_self_merged_mean_coverage,
            "percent_reference_covered" : select_first([percent_reference_covered, 0.0]),

            "intermediate_gapfill_fasta" : scaffold.intermediate_gapfill_fasta,
            "assembly_preimpute_length_unambiguous" : scaffold.assembly_preimpute_length_unambiguous,

            "replicate_concordant_sites" : refine.replicate_concordant_sites,
            "replicate_discordant_snps" : refine.replicate_discordant_snps,
            "replicate_discordant_indels" : refine.replicate_discordant_indels,
            "replicate_discordant_vcf" : refine.replicate_discordant_vcf,

            "isnvsFile" : refine.align_to_self_isnvs_vcf,
            "aligned_bam" : refine.align_to_self_merged_aligned_only_bam,
            "coverage_tsv" : refine.align_to_self_merged_coverage_tsv,
            "read_pairs_aligned" : refine.align_to_self_merged_read_pairs_aligned,
            "bases_aligned" : refine.align_to_self_merged_bases_aligned
        }


        scatter(h in assembly_header) {
            String stat_by_taxon = stats_by_taxon[h]
        }
    }
 
    ### summary stats
    call utils.concatenate {
      input:
        infiles     = [write_tsv([assembly_header]), write_tsv(stat_by_taxon)],
        output_name = "assembly_metadata-~{sample_id}.tsv"
    }

    output {
        Array[Map[String,String]] assembly_stats_by_taxon = stats_by_taxon
        File   assembly_stats_by_taxon_tsv = concatenate.combined

        Int    num_read_groups                       = refine.num_read_groups[0]
        Int    num_libraries                         = refine.num_libraries[0]

        String assembly_method = "viral-ngs/scaffold_and_refine_multitaxa"
        String scaffold_viral_assemble_version       = scaffold.viralngs_version[0]
        String refine_viral_assemble_version         = refine.viral_assemble_version[0]
    }
}
