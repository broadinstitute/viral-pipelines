version 1.0

import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_utils.wdl" as utils

workflow genbank_single {

    meta {
        description: "Prepare assemblies for Genbank submission. This includes annotation by simple coordinate transfer from Genbank annotations and a multiple alignment. See https://viral-pipelines.readthedocs.io/en/latest/ncbi_submission.html for details."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File          assembly_fasta
        String        ref_accessions_colon_delim

        String        biosample_accession
        Int           tax_id
        String        organism_name

        String        email_address # required for fetching data from NCBI APIs

        String?       biosample_attributes_json # if this is used, we will use this first
        File?         biosample_attributes_tsv # if no json, we will read this tsv
        # if both are unspecified, we will fetch from NCBI via biosample_accession
    }

    parameter_meta {
        assembly_fasta: {
          description: "Genome to prepare for Genbank submission. All segments/chromosomes included in one file. Must contain exactly the same number of sequences as reference_accessions.",
          patterns: ["*.fasta"]
        }
        ref_accessions_colon_delim: {
          description: "Reference genome Genbank accessions, each segment/chromosome in the exact same count and order as the segments/chromosomes described in assemblies_fasta. List of accessions should be colon delimited.",
          patterns: ["*.fasta"]
        }
        biosample_attributes_tsv: {
          description: "A post-submission attributes file from NCBI BioSample, which is available at https://submit.ncbi.nlm.nih.gov/subs/ and clicking on 'Download attributes file with BioSample accessions'.",
          patterns: ["*.txt", "*.tsv"]
        }
    }

    # fetch biosample metadata from NCBI if it's not given to us in tsv form
    if(defined(biosample_attributes_json)) {
        call utils.json_dict_to_tsv as biosample_json_to_tsv {
            input:
                json_data = select_first([biosample_attributes_json]),
                out_basename = "biosample_attributes-~{biosample_accession}"
        }
    }
    if(!defined(biosample_attributes_tsv) && !defined(biosample_attributes_json)) {
        call ncbi_tools.fetch_biosamples {
            input:
                biosample_ids = [biosample_accession],
                out_basename = "biosample_attributes-~{biosample_accession}"
        }
    }
    File biosample_attributes = select_first([biosample_json_to_tsv.tsv, biosample_attributes_tsv, fetch_biosamples.biosample_attributes_tsv])

    # Is this a special virus that NCBI handles differently?
    call ncbi.genbank_special_taxa {
      input:
        taxid = tax_id
    }

    # Annotate genes
    ## fetch reference genome sequences and annoations
    call utils.string_split {
        input:
            joined_string = ref_accessions_colon_delim,
            delimiter = ":"
    }
    scatter(segment_acc in string_split.tokens) {
      ## scatter these calls in order to preserve original order
      call ncbi.download_annotations {
        input:
          accessions = [segment_acc],
          emailAddress = email_address,
          combined_out_prefix = segment_acc
      }
    }
    # create genbank source modifier table from biosample metadata
    call ncbi.biosample_to_genbank {
        input:
            biosample_attributes   = biosample_attributes,
            num_segments           = length(string_split.tokens),
            taxid                  = tax_id,
            organism_name_override = organism_name,
            filter_to_accession    = biosample_accession
    }

    ## annotate genes, either by VADR or by naive coordinate liftover
    if(!genbank_special_taxa.vadr_supported) {
      call ncbi.align_and_annot_transfer_single as annot {
        input:
            genome_fasta             = assembly_fasta,
            reference_fastas         = flatten(download_annotations.genomes_fasta),
            reference_feature_tables = flatten(download_annotations.features_tbl)
      }
    }
    if(genbank_special_taxa.vadr_supported) {
      call ncbi.vadr {
        input:
          genome_fasta    = assembly_fasta,
          maxlen          = genbank_special_taxa.max_genome_length,
          vadr_opts       = genbank_special_taxa.vadr_cli_options,
          vadr_model_tar  = genbank_special_taxa.vadr_model_tar,
          mem_size        = genbank_special_taxa.vadr_min_ram_gb
      }
    }
    File feature_tbl   = select_first([vadr.feature_tbl, annot.feature_tbl])

    call ncbi.structured_comments_from_aligned_bam

    if(genbank_special_taxa.table2asn_allowed) {
      call ncbi.table2asn {
        input:
            assembly_fasta          = assembly_fasta,
            annotations_tbl         = feature_tbl,
            source_modifier_table   = biosample_to_genbank.genbank_source_modifier_table,
            structured_comment_file = structured_comments_from_aligned_bam.structured_comment_file,
            organism                = organism_name
      }
    }

    output {
        File          genbank_comment_file   = structured_comments_from_aligned_bam.structured_comment_file
        File          genbank_source_table   = biosample_to_genbank.genbank_source_modifier_table
        File          annotation_tbl         = feature_tbl

        Boolean?      vadr_pass              = vadr.pass
        Array[String] vadr_alerts            = select_first([vadr.alerts, []])

        File?         genbank_submission_sqn = table2asn.genbank_submission_sqn
        File?         genbank_preview_file   = table2asn.genbank_preview_file
        File?         table2asn_val_file     = table2asn.genbank_validation_file
        Array[String] table2asn_errors       = select_first([table2asn.table2asn_errors, []])

        # TO DO: add some conditionally present fields such as File? passing_sqn (only if no errors from vadr or table2asn) or File? passing_non_sqn_zips (for the four special non-table2asn taxa but only if their vadr is passing; zipping up the fsa, cmt, and src files)
    }

}
