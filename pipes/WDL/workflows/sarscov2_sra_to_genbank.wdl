version 1.0

import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_reports.wdl" as reports

import "assemble_refbased.wdl"
import "sarscov2_lineages.wdl"

workflow sarscov2_sra_to_genbank {
    meta {
        description: "Full SARS-CoV-2 analysis workflow starting from SRA data and metadata and performing assembly, spike-in analysis, qc, lineage assignment, and packaging assemblies for data release."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
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

        Int           min_genome_bases = 20000
        Int           min_reads_per_bam = 100
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
            File reads_ubam = Fetch_SRA_to_BAM.reads_ubam
            File biosample_attributes_json = Fetch_SRA_to_BAM.biosample_attributes_json
            String library_strategy = Fetch_SRA_to_BAM.library_strategy
            String biosample_accession = Fetch_SRA_to_BAM.biosample_accession
            Int num_reads = Fetch_SRA_to_BAM.num_reads
            call reports.fastqc {
                input:
                    reads_bam = Fetch_SRA_to_BAM.reads_ubam
            }
            call reports.align_and_count as spikein {
                input:
                    reads_bam = Fetch_SRA_to_BAM.reads_ubam,
                    ref_db = spikein_db
            }
        }
    }

    ### gather data by biosample
    call ncbi_tools.group_sra_bams_by_biosample {
        input:
            biosamples                 = select_all(biosample_accession),
            bam_filepaths              = select_all(reads_ubam),
            biosample_attributes_jsons = select_all(biosample_attributes_json),
            library_strategies         = select_all(library_strategy)
    }

    ### assembly and analyses per biosample
    scatter(samn_bam in zip(group_sra_bams_by_biosample.biosample_accessions, group_sra_bams_by_biosample.grouped_bam_filepaths)) {
        Boolean ampseq = (group_sra_bams_by_biosample.samn_to_library_strategy[samn_bam.left] == "AMPLICON")
        String orig_name = group_sra_bams_by_biosample.samn_to_attributes[samn_bam.left]["sample_name"]

        # assemble genome
        if (ampseq) {
            String trim_coords_bed = amplicon_bed_default
        }
        call assemble_refbased.assemble_refbased {
            input:
                reads_unmapped_bams = samn_bam.right,
                reference_fasta = reference_fasta,
                sample_name = samn_bam.left,
                aligner = "minimap2",
                skip_mark_dupes = ampseq,
                trim_coords_bed = trim_coords_bed,
                min_coverage = if ampseq then 20 else 3
        }

        # for genomes that somewhat assemble
        if (assemble_refbased.assembly_length_unambiguous >= min_genome_bases) {
            call ncbi.rename_fasta_header {
              input:
                genome_fasta = assemble_refbased.assembly_fasta,
                new_name = orig_name
            }

            File passing_assemblies = rename_fasta_header.renamed_fasta
            String passing_assembly_ids = orig_name
            Array[String] assembly_cmt = [orig_name, "Broad viral-ngs v. " + assemble_refbased.align_to_ref_viral_core_version, assemble_refbased.assembly_mean_coverage]

            # lineage assignment
            call sarscov2_lineages.sarscov2_lineages {
                input:
                    genome_fasta = passing_assemblies
            }

            # VADR annotation & QC
            call ncbi.vadr {
              input:
                genome_fasta = passing_assemblies
            }
            if (vadr.num_alerts==0) {
              File submittable_genomes = passing_assemblies
              String submittable_id = orig_name
            }
            if (vadr.num_alerts>0) {
              String failed_annotation_id = orig_name
            }
        }
        if (assemble_refbased.assembly_length_unambiguous < min_genome_bases) {
            String failed_assembly_id = orig_name
        }

        Map[String,String?] assembly_stats = {
            'sample_orig': orig_name,
            'sample': samn_bam.left,
            'assembly_mean_coverage': assemble_refbased.assembly_mean_coverage,
            'nextclade_clade':   sarscov2_lineages.nextclade_clade,
            'nextclade_aa_subs': sarscov2_lineages.nextclade_aa_subs,
            'nextclade_aa_dels': sarscov2_lineages.nextclade_aa_dels,
            'pango_lineage':     sarscov2_lineages.pango_lineage
        }
        Map[String,File?] assembly_files = {
            'assembly_fasta':           assemble_refbased.assembly_fasta,
            'coverage_plot':            assemble_refbased.align_to_ref_merged_coverage_plot,
            'aligned_bam':              assemble_refbased.align_to_ref_merged_aligned_trimmed_only_bam,
            'replicate_discordant_vcf': assemble_refbased.replicate_discordant_vcf,
            'nextclade_tsv': sarscov2_lineages.nextclade_tsv,
            'pangolin_csv':  sarscov2_lineages.pangolin_csv,
            'vadr_tgz': vadr.outputs_tgz
        }
        Map[String,Int?] assembly_metrics = {
            'assembly_length_unambiguous': assemble_refbased.assembly_length_unambiguous,
            'dist_to_ref_snps':            assemble_refbased.dist_to_ref_snps,
            'dist_to_ref_indels':          assemble_refbased.dist_to_ref_indels,
            'replicate_concordant_sites':  assemble_refbased.replicate_concordant_sites,
            'replicate_discordant_snps':   assemble_refbased.replicate_discordant_snps,
            'replicate_discordant_indels': assemble_refbased.replicate_discordant_indels,
            'num_read_groups':             assemble_refbased.num_read_groups,
            'num_libraries':               assemble_refbased.num_libraries,
            'vadr_num_alerts': vadr.num_alerts
        }

    }

    ### prep genbank submission
    call nextstrain.concatenate as submit_genomes {
      input:
        infiles = select_all(submittable_genomes),
        output_name = "assemblies.fasta"
    }
    call ncbi.biosample_to_genbank {
      input:
        biosample_attributes = group_sra_bams_by_biosample.biosample_attributes_tsv,
        num_segments = 1,
        taxid = taxid,
        filter_to_ids = write_lines(select_all(submittable_id))
    }
    call ncbi.structured_comments {
      input:
        assembly_stats_tsv = write_tsv(flatten([[['SeqID','Assembly Method','Coverage']],select_all(assembly_cmt)])),
        filter_to_ids = write_lines(select_all(submittable_id))
    }
    call ncbi.package_genbank_ftp_submission {
      input:
        sequences_fasta = submit_genomes.combined,
        source_modifier_table = biosample_to_genbank.genbank_source_modifier_table,
        structured_comment_table = structured_comments.structured_comment_table
    }

    ### prep gisaid submission
    call ncbi.prefix_fasta_header as prefix_gisaid {
      input:
        genome_fasta = submit_genomes.combined,
        prefix = gisaid_prefix
    }
    call ncbi.gisaid_meta_prep {
      input:
        source_modifier_table = biosample_to_genbank.genbank_source_modifier_table,
        structured_comments = structured_comments.structured_comment_table,
        out_name = "gisaid_meta.tsv",
        strict = false
    }

    #### summary stats
    call reports.MultiQC as multiqc_raw {
        input:
            input_files = select_all(fastqc.fastqc_zip),
            file_name   = "multiqc-raw.html"
    }
    call reports.align_and_count_summary as spike_summary {
        input:
            counts_txt = select_all(spikein.report)
    }

    output {
        Array[File] raw_reads_unaligned_bams     = select_all(reads_ubam)
        Array[File] cleaned_reads_unaligned_bams = select_all(reads_ubam)

        File        biosample_attributes = group_sra_bams_by_biosample.biosample_attributes_tsv

        Array[Int]  read_counts_raw       = select_all(num_reads)
        Array[Int]  read_counts_depleted  = select_all(num_reads)

        Array[File] assemblies_fasta = assemble_refbased.assembly_fasta
        Array[File] passing_assemblies_fasta = select_all(passing_assemblies)
        Array[File] submittable_assemblies_fasta = select_all(submittable_genomes)

        File        multiqc_report_raw     = multiqc_raw.multiqc_report
        File        spikein_counts         = spike_summary.count_summary

        Array[Map[String,String?]] per_assembly_stats = assembly_stats
        Array[Map[String,File?]]   per_assembly_files = assembly_files
        Array[Map[String,Int?]]    per_assembly_metrics = assembly_metrics

        File submission_zip = package_genbank_ftp_submission.submission_zip
        File submission_xml = package_genbank_ftp_submission.submission_xml
        File submit_ready   = package_genbank_ftp_submission.submit_ready
        Array[File] vadr_outputs = select_all(vadr.outputs_tgz)
        File genbank_source_table = biosample_to_genbank.genbank_source_modifier_table

        File gisaid_fasta = prefix_gisaid.renamed_fasta
        File gisaid_meta_tsv = gisaid_meta_prep.meta_tsv

        Array[String] assembled_ids = select_all(passing_assembly_ids)
        Array[String] submittable_ids = select_all(submittable_id)
        Array[String] failed_assembly_ids = select_all(failed_assembly_id)
        Array[String] failed_annotation_ids = select_all(failed_annotation_id)
        Int           num_read_files = length(Fetch_SRA_to_BAM.reads_ubam)
        Int           num_assembled = length(select_all(passing_assemblies))
        Int           num_failed_assembly = length(select_all(failed_assembly_id))
        Int           num_submittable = length(select_all(submittable_id))
        Int           num_failed_annotation = length(select_all(failed_annotation_id))
        Int           num_samples = length(group_sra_bams_by_biosample.biosample_accessions)
    }
}
