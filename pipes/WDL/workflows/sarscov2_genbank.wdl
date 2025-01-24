version 1.0

import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_utils.wdl" as utils

workflow sarscov2_genbank {

    meta {
        description: "Prepare SARS-CoV-2 assemblies for Genbank submission. This includes QC checks with NCBI's VADR tool and filters out genomes that do not pass its tests."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+  assemblies_fasta

        String?       author_list # of the form "Lastname,A.B., Lastname,C.,"; optional alternative to names in author_sbt_defaults_yaml
        File          author_sbt_defaults_yaml # defaults to fill in for author_sbt file (including both author and non-author fields)
        File          author_sbt_j2_template
        File          biosample_attributes
        File          assembly_stats_tsv
        File?         fasta_rename_map

        Int           min_genome_bases = 15000
        Int           max_vadr_alerts = 0

        Int           taxid = 2697049
        String        gisaid_prefix = 'hCoV-19/'
    }

    parameter_meta {
        assemblies_fasta: {
          description: "Genomes to prepare for Genbank submission. One file per genome: all segments/chromosomes included in one file. All fasta files must contain exactly the same number of sequences as reference_fasta (which must equal the number of files in reference_annot_tbl).",
          patterns: ["*.fasta"]
        }
        author_list: {
          description: "A string containing a space-delimited list with of author surnames separated by first name and (optional) middle initial. Ex. 'Lastname,Firstname, Last-hypenated,First,M., Last,F.'"
        }
        author_sbt_defaults_yaml: {
          description: "A YAML file with default values to use for the submitter, submitter affiliation, and author affiliation. Optionally including authors at the start and end of the author_list. Example: gs://pathogen-public-dbs/other-related/default_sbt_values.yaml",
          patterns: ["*.yaml","*.yml"]
        }
        author_sbt_j2_template: {
          description: "A jinja2-format template for the sbt file expected by NCBI. Example: gs://pathogen-public-dbs/other-related/author_template.sbt.j2"
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
          String fasta_basename = basename(assembly, ".fasta")
          call ncbi.rename_fasta_header {
            input:
              genome_fasta = assembly,
              new_name     = read_map(select_first([fasta_rename_map]))[fasta_basename]
          }
        }
        call reports.assembly_bases {
          input:
            fasta = assembly
        }
        File renamed_assembly = select_first([rename_fasta_header.renamed_fasta, assembly])
        call ncbi.vadr {
          input:
            genome_fasta = renamed_assembly,
            vadr_opts = "--glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn",
            minlen = 50,
            maxlen = 30000
        }
        if (assembly_bases.assembly_length_unambiguous >= min_genome_bases) {
          if (vadr.num_alerts <= max_vadr_alerts) {
            File passing_assemblies = renamed_assembly
          }
          if (vadr.num_alerts > max_vadr_alerts) {
            File weird_assemblies = renamed_assembly
          }
       }
    }

    # prep the good ones
    call utils.concatenate as passing_fasta {
      input:
        infiles     = select_all(passing_assemblies),
        output_name = "assemblies-passing.fasta"
    }
    call utils.fasta_to_ids as passing_ids {
      input:
        sequences_fasta = passing_fasta.combined
    }

    # prep the weird ones
    call utils.concatenate as weird_fasta {
      input:
        infiles     = select_all(weird_assemblies),
        output_name = "assemblies-weird.fasta"
    }
    call utils.fasta_to_ids as weird_ids {
      input:
        sequences_fasta = weird_fasta.combined
    }

    # package genbank
    call ncbi.biosample_to_genbank as passing_source_modifiers {
      input:
        biosample_attributes = biosample_attributes,
        num_segments         = 1,
        taxid                = taxid,
        filter_to_ids        = passing_ids.ids_txt
    }
    call ncbi.structured_comments as passing_structured_cmt {
      input:
        assembly_stats_tsv = assembly_stats_tsv,
        filter_to_ids      = passing_ids.ids_txt
    }
    call ncbi.generate_author_sbt_file as generate_author_sbt {
        input:
            author_list   = author_list,
            defaults_yaml = author_sbt_defaults_yaml,
            j2_template   = author_sbt_j2_template
    }
    call ncbi.package_sc2_genbank_ftp_submission as passing_package_genbank {
      input:
        sequences_fasta          = passing_fasta.combined,
        source_modifier_table    = passing_source_modifiers.genbank_source_modifier_table,
        author_template_sbt      = generate_author_sbt.sbt_file,
        structured_comment_table = passing_structured_cmt.structured_comment_table
    }

    # translate to gisaid
    call ncbi.prefix_fasta_header as passing_prefix_gisaid {
      input:
        genome_fasta = passing_fasta.combined,
        prefix       = gisaid_prefix,
        out_basename = "gisaid-passing-sequences"
    }
    call ncbi.gisaid_meta_prep as passing_gisaid_meta {
      input:
        source_modifier_table = passing_source_modifiers.genbank_source_modifier_table,
        structured_comments   = passing_structured_cmt.structured_comment_table,
        fasta_filename        = "gisaid-passing-sequences.fasta",
        out_name              = "gisaid-passing-meta.csv"
    }

    # package genbank
    call ncbi.biosample_to_genbank as weird_source_modifiers {
      input:
        biosample_attributes = biosample_attributes,
        num_segments         = 1,
        taxid                = taxid,
        filter_to_ids        = weird_ids.ids_txt
    }
    call ncbi.structured_comments as weird_structured_cmt {
      input:
        assembly_stats_tsv = assembly_stats_tsv,
        filter_to_ids      = weird_ids.ids_txt
    }
    call ncbi.package_sc2_genbank_ftp_submission as weird_package_genbank {
      input:
        sequences_fasta          = weird_fasta.combined,
        source_modifier_table    = weird_source_modifiers.genbank_source_modifier_table,
        author_template_sbt      = generate_author_sbt.sbt_file,
        structured_comment_table = weird_structured_cmt.structured_comment_table
    }

    # translate to gisaid
    call ncbi.prefix_fasta_header as weird_prefix_gisaid {
      input:
        genome_fasta = weird_fasta.combined,
        prefix       = gisaid_prefix,
        out_basename = "gisaid-weird-sequences"
    }
    call ncbi.gisaid_meta_prep as weird_gisaid_meta {
      input:
        source_modifier_table = weird_source_modifiers.genbank_source_modifier_table,
        structured_comments   = weird_structured_cmt.structured_comment_table,
        fasta_filename        = "gisaid-weird-sequences.fasta",
        out_name              = "gisaid-weird-meta.csv"
    }

    output {
        File        submission_zip        = passing_package_genbank.submission_zip
        File        submission_xml        = passing_package_genbank.submission_xml
        File        submit_ready          = passing_package_genbank.submit_ready
        
        Int         num_successful        = length(select_all(passing_assemblies))
        Int         num_weird             = length(select_all(weird_assemblies))
        Int         num_input             = length(assemblies_fasta)
        
        Array[File] vadr_outputs          = vadr.outputs_tgz
        
        File        gisaid_fasta          = passing_prefix_gisaid.renamed_fasta
        File        gisaid_meta_csv       = passing_gisaid_meta.meta_csv
        
        File        weird_genbank_zip     = weird_package_genbank.submission_zip
        File        weird_genbank_xml     = weird_package_genbank.submission_xml
        File        weird_gisaid_fasta    = weird_prefix_gisaid.renamed_fasta
        File        weird_gisaid_meta_csv = weird_gisaid_meta.meta_csv
    }

}
