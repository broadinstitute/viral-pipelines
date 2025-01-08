version 1.0

import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_reports.wdl" as reports
import "../tasks/tasks_utils.wdl" as utils

workflow genbank_single {

    meta {
        description: "Prepare assemblies for Genbank submission. This includes annotation by simple coordinate transfer from Genbank annotations and a multiple alignment. See https://viral-pipelines.readthedocs.io/en/latest/ncbi_submission.html for details."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        File            assembly_fasta
        Array[String]+  reference_accessions ### TO DO: take a colon-delim list instead?

        File          aligned_reads_bam
        String        biosample_accession
        Int           tax_id
        String        organism_name
        String        assembly_method
        String        assembly_method_version

        String        email_address # required for fetching data from NCBI APIs
        String?       author_list # of the form "Lastname,A.B., Lastname,C.,"; optional alternative to names in author_sbt_defaults_yaml
        File?         author_sbt_defaults_yaml # defaults to fill in for author_sbt file (including both author and non-author fields)
        File          author_sbt_j2_template # an sbt file (optionally) with Jinja2 variables filled in based on author_sbt_defaults_yaml if provided
        File?         biosample_attributes_tsv # if empty, we will fetch from NCBI via accession
        String?       comment
        String        molType='cRNA'
    }

    parameter_meta {
        assembly_fasta: {
          description: "Genome to prepare for Genbank submission. All segments/chromosomes included in one file. Must contain exactly the same number of sequences as reference_accessions.",
          patterns: ["*.fasta"]
        }
        reference_accessions: {
          description: "Reference genome Genbank accessions, each segment/chromosome in the exact same count and order as the segments/chromosomes described in assemblies_fasta.",
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
          description: "an sbt file (optionally) with Jinja2 variables to be filled in based on values present in author_sbt_defaults_yaml, if provided. If author_list is blank and author_sbt_defaults_yaml is not provided (or is blank), this file is passed through verbatim. Example: gs://pathogen-public-dbs/other-related/author_template.sbt.j2"
        }
        biosample_attributes_tsv: {
          description: "A post-submission attributes file from NCBI BioSample, which is available at https://submit.ncbi.nlm.nih.gov/subs/ and clicking on 'Download attributes file with BioSample accessions'.",
          patterns: ["*.txt", "*.tsv"]
        }
        coverage_table: {
          description: "A two column tab text file mapping sample IDs (first column) to average sequencing coverage (second column, floating point number).",
          patterns: ["*.txt", "*.tsv"],
          category: "common"
        }
        molType: {
          description: "The type of molecule being described. This defaults to 'cRNA' as this pipeline is most commonly used for viral submissions, but any value allowed by the INSDC controlled vocabulary may be used here. Valid values are described at http://www.insdc.org/controlled-vocabulary-moltype-qualifier",
          category: "common"
        }
        comment: {
          description: "Optional comments that can be displayed in the COMMENT section of the Genbank record. This may include any disclaimers about assembly quality or notes about pre-publication availability or requests to discuss pre-publication use with authors."
        }

    }

    # fetch biosample metadata from NCBI if it's not given to us in tsv form
    if(!defined(biosample_attributes_tsv)) {
        call ncbi_tools.fetch_biosamples {
            input:
                biosample_ids = [biosample_accession],
                out_basename = "biosample_attributes-~{biosample_accession}"
        }
    }
    File biosample_attributes = select_first([biosample_attributes_tsv, fetch_biosamples.biosample_attributes_tsv])

    scatter(segment_acc in reference_accessions) {
      # scatter these calls in order to preserve original order
      call ncbi.download_annotations {
        input:
          accessions = [segment_acc],
          emailAddress = email_address,
          combined_out_prefix = segment_acc
      }
    }

    # extract info from aligned bams (coverage, seq tech)
    call ncbi.sequencing_platform_from_bam {
      input:
        bam = aligned_reads_bam
    }
    call reports.coverage_report {
      input:
        mapped_bams = [aligned_reads_bam],
        mapped_bam_idx = []
    }
    call utils.tsv_drop_cols as coverage_two_col {
      input:
        in_tsv = coverage_report.coverage_report,
        drop_cols = ["aln2self_cov_median", "aln2self_cov_mean_non0", "aln2self_cov_1X", "aln2self_cov_5X", "aln2self_cov_20X", "aln2self_cov_100X"]
    }

    # create genbank source modifier table from biosample metadata
    call ncbi.biosample_to_genbank {
        input:
            biosample_attributes = biosample_attributes,
            num_segments         = length(reference_accessions),
            taxid                = tax_id
    }

    # TO DO: for VADR-based viruses, call VADR instead of align_and_annot_transfer_single
    # RSV A & B, MPXV plus also the special ones: SC2, Flu A/B/C, Noro, DENV
    call ncbi.align_and_annot_transfer_single as annot {
        input:
            genome_fasta             = assembly_fasta,
            reference_fastas         = flatten(download_annotations.genomes_fasta),
            reference_feature_tables = flatten(download_annotations.features_tbl)
    }

    call ncbi.generate_author_sbt_file as generate_author_sbt {
        input:
            author_list   = author_list,
            defaults_yaml = author_sbt_defaults_yaml,
            j2_template   = author_sbt_j2_template
    }

    # TO DO: create entirely different branch of pipeline for special viruses (SC2, Flu A/B/C, Noro, DENV)
    # those viruses need only: fasta, source modifier table, "source info"?, "references"?, "sequence processing"?
    # concatenated fastas are useful for bulk submission
    # table2asn outputs are not accepted by Genbank for these
    call ncbi.prepare_genbank_single as prep_genbank {
        input:
            assemblies_fasta   = [assembly_fasta],
            annotations_tbl    = annot.genome_per_chr_tbls,
            authors_sbt        = generate_author_sbt.sbt_file,
            biosampleMap       = biosample_to_genbank.biosample_map,
            genbankSourceTable = biosample_to_genbank.genbank_source_modifier_table,
            coverage_table     = coverage_two_col.out_tsv,
            sequencingTech     = sequencing_platform_from_bam.genbank_sequencing_technology,
            comment            = comment,
            organism           = organism_name,
            molType            = molType,
            assembly_method    = assembly_method,
            assembly_method_version = assembly_method_version
    }

    output {
        File        submission_zip         = prep_genbank.submission_zip
        File        archive_zip            = prep_genbank.archive_zip
        File        errorSummary           = prep_genbank.errorSummary
        
        File        genbank_source_table   = biosample_to_genbank.genbank_source_modifier_table
        
        Array[File] transferred_annot_tbls = annot.genome_per_chr_tbls
        Array[File] genbank_preview_files  = prep_genbank.genbank_preview_files
        Array[File] validation_files       = prep_genbank.validation_files
        
        String      viral_phylo_version    = prep_genbank.viralngs_version
    }

}
