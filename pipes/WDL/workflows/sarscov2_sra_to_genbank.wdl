version 1.1

#DX_SKIP_WORKFLOW

import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_utils.wdl" as utils

import "assemble_refbased.wdl"
import "sarscov2_lineages.wdl"

workflow sarscov2_sra_to_genbank {
    meta {
        description: "Full SARS-CoV-2 analysis workflow starting from SRA data and metadata and performing assembly, spike-in analysis, qc, lineage assignment, and packaging assemblies for data release."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    parameter_meta {
        SRA_accessions: {
            description: "SRA Run accession numbers (SRR####)."
        }
        reference_fasta: {
            description: "Reference genome to align reads to.",
            patterns: ["*.fasta"]
        }
        amplicon_bed_default: {
            description: "Amplicon primers to trim in reference coordinate space (0-based BED format). Will only be used if the SRA Experiment Design has a Library Strategy set to AMPLICON (will be ignored on all other values).",
            patterns: ["*.bed"]
        }
    }

    input {
        Array[String] SRA_accessions

        File          reference_fasta
        File          amplicon_bed_default
        File          spikein_db

        Int           min_genome_bases = 15000
        Int           min_reads_per_bam = 100
        Int           max_vadr_alerts = 0
    }
    Int     taxid = 2697049
    String  gisaid_prefix = 'hCoV-19/'

    ### retrieve reads and annotation from SRA
    scatter(sra_id in SRA_accessions) {
        call ncbi_tools.Fetch_SRA_to_BAM {
            input:
                SRA_ID = sra_id
        }
        if (Fetch_SRA_to_BAM.num_reads >= min_reads_per_bam) {
            File reads_ubam                = Fetch_SRA_to_BAM.reads_ubam
            File biosample_attributes_json = Fetch_SRA_to_BAM.biosample_attributes_json
            String library_strategy        = Fetch_SRA_to_BAM.library_strategy
            String biosample_accession     = Fetch_SRA_to_BAM.biosample_accession
            Int num_reads                  = Fetch_SRA_to_BAM.num_reads
            String seq_platform            = Fetch_SRA_to_BAM.sequencing_platform
            call reports.align_and_count as spikein {
                input:
                    reads_bam = Fetch_SRA_to_BAM.reads_ubam,
                    ref_db    = spikein_db
            }
        }
    }

    ### gather data by biosample
    call ncbi_tools.group_sra_bams_by_biosample {
        input:
            biosamples                 = select_all(biosample_accession),
            bam_filepaths              = select_all(reads_ubam),
            biosample_attributes_jsons = select_all(biosample_attributes_json),
            library_strategies         = select_all(library_strategy),
            seq_platforms              = select_all(seq_platform)
    }

    ### assembly and analyses per biosample
    scatter(samn_bam in zip(group_sra_bams_by_biosample.biosample_accessions, group_sra_bams_by_biosample.grouped_bam_filepaths)) {
        Boolean ampseq   = (group_sra_bams_by_biosample.samn_to_library_strategy[samn_bam.left] == "AMPLICON")
        String orig_name = group_sra_bams_by_biosample.samn_to_attributes[samn_bam.left]["sample_name"]

        # assemble genome
        if (ampseq) {
            String trim_coords_bed = amplicon_bed_default
        }
        call assemble_refbased.assemble_refbased {
            input:
                reads_unmapped_bams = samn_bam.right,
                reference_fasta     = reference_fasta,
                sample_name         = samn_bam.left,
                aligner             = "minimap2",
                skip_mark_dupes     = ampseq,
                trim_coords_bed     = trim_coords_bed,
                min_coverage        = if ampseq then 20 else 3
        }

        # for genomes that somewhat assemble
        if (assemble_refbased.assembly_length_unambiguous >= min_genome_bases) {
            call ncbi.rename_fasta_header {
              input:
                genome_fasta = assemble_refbased.assembly_fasta,
                new_name     = orig_name
            }

            File passing_assemblies     = rename_fasta_header.renamed_fasta
            String passing_assembly_ids = orig_name
            Array[String] assembly_cmt  = [orig_name, "Broad viral-ngs v. " + assemble_refbased.align_to_ref_viral_core_version, assemble_refbased.assembly_mean_coverage, group_sra_bams_by_biosample.samn_to_sequencing_platform[samn_bam.left]]

            # lineage assignment
            call sarscov2_lineages.sarscov2_lineages {
                input:
                    genome_fasta = passing_assemblies
            }

            # VADR annotation & QC
            call ncbi.vadr {
              input:
                genome_fasta = passing_assemblies,
                vadr_opts = "--glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn",
                minlen = 50,
                maxlen = 30000
            }
            if (vadr.num_alerts<=max_vadr_alerts) {
              File submittable_genomes = passing_assemblies
              String submittable_id    = orig_name
            }
            if (vadr.num_alerts>max_vadr_alerts) {
              String failed_annotation_id = orig_name
            }
        }
        if (assemble_refbased.assembly_length_unambiguous < min_genome_bases) {
            String failed_assembly_id = orig_name
        }

        Array[String] assembly_tsv_row = [
            orig_name,
            samn_bam.left,
            basename(select_first([trim_coords_bed, ""])),
            assemble_refbased.assembly_mean_coverage,
            assemble_refbased.assembly_length_unambiguous,
            select_first([sarscov2_lineages.nextclade_clade, ""]),
            select_first([sarscov2_lineages.nextclade_aa_subs, ""]),
            select_first([sarscov2_lineages.nextclade_aa_dels, ""]),
            select_first([sarscov2_lineages.pango_lineage, ""]),
            assemble_refbased.dist_to_ref_snps,
            assemble_refbased.dist_to_ref_indels,
            select_first([vadr.num_alerts, ""]),
            assemble_refbased.assembly_fasta,
            assemble_refbased.align_to_ref_merged_coverage_plot,
            assemble_refbased.align_to_ref_merged_aligned_trimmed_only_bam,
            assemble_refbased.replicate_discordant_vcf,
            select_first([sarscov2_lineages.nextclade_tsv, ""]),
            select_first([sarscov2_lineages.nextclade_json, ""]),
            select_first([sarscov2_lineages.pango_lineage_report, ""]),
            select_first([vadr.outputs_tgz, ""]),
            assemble_refbased.replicate_concordant_sites,
            assemble_refbased.replicate_discordant_snps,
            assemble_refbased.replicate_discordant_indels,
            assemble_refbased.num_read_groups,
            assemble_refbased.num_libraries,
            assemble_refbased.align_to_ref_merged_reads_aligned,
            assemble_refbased.align_to_ref_merged_bases_aligned,
        ]
    }

    Array[String] assembly_tsv_header = [
        'sample', 'sample_sanitized', 'amplicon_set', 'assembly_mean_coverage', 'assembly_length_unambiguous',
        'nextclade_clade', 'nextclade_aa_subs', 'nextclade_aa_dels', 'pango_lineage',
        'dist_to_ref_snps', 'dist_to_ref_indels', 'vadr_num_alerts',
        'assembly_fasta', 'coverage_plot', 'aligned_bam', 'replicate_discordant_vcf',
        'nextclade_tsv', 'nextclade_json', 'pangolin_csv', 'vadr_tgz',
        'replicate_concordant_sites', 'replicate_discordant_snps', 'replicate_discordant_indels', 'num_read_groups', 'num_libraries',
        'align_to_ref_merged_reads_aligned', 'align_to_ref_merged_bases_aligned',
        ]

    call utils.concatenate as assembly_meta_tsv {
      input:
        infiles = [write_tsv([assembly_tsv_header]), write_tsv(assembly_tsv_row)],
        output_name = "assembly_metadata.tsv"
    }

    ### prep genbank submission
    call utils.concatenate as submit_genomes {
      input:
        infiles     = select_all(submittable_genomes),
        output_name = "assemblies.fasta"
    }
    call ncbi.biosample_to_genbank {
      input:
        biosample_attributes = group_sra_bams_by_biosample.biosample_attributes_tsv,
        num_segments         = 1,
        taxid                = taxid,
        filter_to_ids        = write_lines(select_all(submittable_id))
    }
    call ncbi.structured_comments {
      input:
        assembly_stats_tsv = write_tsv(flatten([[['SeqID','Assembly Method','Coverage','Sequencing Technology']],select_all(assembly_cmt)])),
        filter_to_ids      = write_lines(select_all(submittable_id))
    }
    call ncbi.package_special_genbank_ftp_submission as package_genbank_ftp_submission {
      input:
        sequences_fasta          = submit_genomes.combined,
        source_modifier_table    = biosample_to_genbank.genbank_source_modifier_table,
        structured_comment_table = structured_comments.structured_comment_table
    }

    ### prep gisaid submission
    call ncbi.prefix_fasta_header as prefix_gisaid {
      input:
        genome_fasta = submit_genomes.combined,
        prefix       = gisaid_prefix
    }
    call ncbi.gisaid_meta_prep {
      input:
        source_modifier_table = biosample_to_genbank.genbank_source_modifier_table,
        structured_comments   = structured_comments.structured_comment_table,
        out_name              = "gisaid_meta.csv",
        strict                = false
    }

    #### summary stats
    call reports.multiqc_from_bams as multiqc_raw {
        input:
            input_bams   = select_all(reads_ubam),
            out_basename = "multiqc-raw"
    }
    call reports.align_and_count_summary as spike_summary {
        input:
            counts_txt = select_all(spikein.report)
    }

    output {
        Array[File]   raw_reads_unaligned_bams     = select_all(reads_ubam)
        Array[File]   cleaned_reads_unaligned_bams = select_all(reads_ubam)
        
        File          biosample_attributes         = group_sra_bams_by_biosample.biosample_attributes_tsv
        
        Array[Int]    read_counts_raw              = select_all(num_reads)
        Array[Int]    read_counts_depleted         = select_all(num_reads)
        
        Array[Int]    primer_trimmed_read_count    = flatten(assemble_refbased.primer_trimmed_read_count)
        Array[Float]  primer_trimmed_read_percent  = flatten(assemble_refbased.primer_trimmed_read_percent)
        
        Array[File]   assemblies_fasta             = assemble_refbased.assembly_fasta
        Array[File]   passing_assemblies_fasta     = select_all(passing_assemblies)
        Array[File]   submittable_assemblies_fasta = select_all(submittable_genomes)
        
        File          multiqc_report_raw           = multiqc_raw.multiqc_report
        File          spikein_counts               = spike_summary.count_summary
        
        File          assembly_stats_tsv           = assembly_meta_tsv.combined
        
        File          submission_zip               = package_genbank_ftp_submission.submission_zip
        File          submission_xml               = package_genbank_ftp_submission.submission_xml
        File          submit_ready                 = package_genbank_ftp_submission.submit_ready
        Array[File]   vadr_outputs                 = select_all(vadr.outputs_tgz)
        File          genbank_source_table         = biosample_to_genbank.genbank_source_modifier_table
        
        File          gisaid_fasta                 = prefix_gisaid.renamed_fasta
        File          gisaid_meta_csv              = gisaid_meta_prep.meta_csv
        
        Array[String] assembled_ids                = select_all(passing_assembly_ids)
        Array[String] submittable_ids              = select_all(submittable_id)
        Array[String] failed_assembly_ids          = select_all(failed_assembly_id)
        Array[String] failed_annotation_ids        = select_all(failed_annotation_id)
        Int           num_read_files               = length(Fetch_SRA_to_BAM.reads_ubam)
        Int           num_assembled                = length(select_all(passing_assemblies))
        Int           num_failed_assembly          = length(select_all(failed_assembly_id))
        Int           num_submittable              = length(select_all(submittable_id))
        Int           num_failed_annotation        = length(select_all(failed_annotation_id))
        Int           num_samples                  = length(group_sra_bams_by_biosample.biosample_accessions)
    }
}
