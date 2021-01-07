version 1.0

import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_nextstrain.wdl" as nextstrain
import "../tasks/tasks_reports.wdl" as reports

workflow sarscov2_genbank {

    meta {
        description: "Prepare SARS-CoV-2 assemblies for Genbank submission. This includes QC checks with NCBI's VADR tool and filters out genomes that do not pass its tests."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+  assemblies_fasta

        File          authors_sbt
        File          biosample_attributes
        File          assembly_stats_tsv
        File?         fasta_rename_map

        Int           taxid = 2697049
        Int           min_genome_bases = 20000
    }

    parameter_meta {
        assemblies_fasta: {
          description: "Genomes to prepare for Genbank submission. One file per genome: all segments/chromosomes included in one file. All fasta files must contain exactly the same number of sequences as reference_fasta (which must equal the number of files in reference_annot_tbl).",
          patterns: ["*.fasta"]
        }
        authors_sbt: {
          description: "A genbank submission template file (SBT) with the author list, created at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/",
          patterns: ["*.sbt"]
        }
        biosample_attributes: {
          description: "A post-submission attributes file from NCBI BioSample, which is available at https://submit.ncbi.nlm.nih.gov/subs/ and clicking on 'Download attributes file with BioSample accessions'.",
          patterns: ["*.txt", "*.tsv"]
        }
        assembly_stats_tsv: {
          description: "A four column tab text file with one row per sequence and the following header columns: SeqID, Assembly Method, Coverage, Sequencing Technology",
          patterns: ["*.txt", "*.tsv"]
        }
    }

    scatter(assembly in assemblies_fasta) {
        if(defined(fasta_rename_map)) {
          call ncbi.lookup_table_by_filename {
            input:
              id = basename(assembly),
              mapping_tsv = select_first([fasta_rename_map])
          }
          call ncbi.rename_fasta {
            input:
              genome_fasta = assembly,
              new_name = lookup_table_by_filename.value
          }
        }
        File renamed_assembly = select_first([rename_fasta.renamed_fasta, assembly])
        call reports.assembly_bases {
          input:
            fasta = renamed_assembly
        }
        call ncbi.vadr {
          input:
            genome_fasta = renamed_assembly
        }
        if ( (vadr.num_alerts==0) && (assembly_bases.assembly_length_unambiguous >= min_genome_bases) ) {
          File passing_assemblies = renamed_assembly
        }
    }

    call nextstrain.concatenate {
      input:
        infiles = select_all(passing_assemblies),
        output_name = "assemblies.fasta"
    }
    call nextstrain.fasta_to_ids {
      input:
        sequences_fasta = concatenate.combined
    }

    call ncbi.biosample_to_genbank {
      input:
        biosample_attributes = biosample_attributes,
        num_segments = 1,
        taxid = taxid,
        filter_to_ids = fasta_to_ids.ids_txt
    }

    call ncbi.structured_comments {
      input:
        assembly_stats_tsv = assembly_stats_tsv,
        filter_to_ids = fasta_to_ids.ids_txt
    }

    call ncbi.package_genbank_ftp_submission {
      input:
        sequences_fasta = concatenate.combined,
        source_modifier_table = biosample_to_genbank.genbank_source_modifier_table,
        author_template_sbt = authors_sbt,
        structured_comment_table = structured_comments.structured_comment_table
    }

    output {
        File submission_zip = package_genbank_ftp_submission.submission_zip
        File submission_xml    = package_genbank_ftp_submission.submission_xml
        File submit_ready   = package_genbank_ftp_submission.submit_ready

        Int  num_successful = length(select_all(passing_assemblies))
        Int  num_input = length(assemblies_fasta)

        Array[File] vadr_outputs = vadr.outputs_tgz

        File biosample_map = biosample_to_genbank.biosample_map
        File genbank_source_table = biosample_to_genbank.genbank_source_modifier_table
    }

}
