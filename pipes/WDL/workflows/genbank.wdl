version 1.0

import "../tasks/tasks_interhost.wdl" as interhost
import "../tasks/tasks_ncbi.wdl" as ncbi

workflow genbank {

    meta {
        description: "Prepare assemblies for Genbank submission. This includes annotation by simple coordinate transfer from Genbank annotations and a multiple alignment. See https://viral-pipelines.readthedocs.io/en/latest/ncbi_submission.html for details."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+  reference_fastas
        Array[File]+  reference_feature_tables
        Array[File]+  assemblies_fasta

        String?       author_list # of the form "Lastname,A.B., Lastname,C.,"; optional alternative to names in author_sbt_defaults_yaml
        File          author_sbt_defaults_yaml # defaults to fill in for author_sbt file (including both author and non-author fields)
        File          author_sbt_j2_template
        File          biosample_attributes
        Int           taxid
        File?         coverage_table
        String?       sequencingTech
        String?       comment
        String?       organism
        String?       molType='cRNA'
    }

    parameter_meta {
        assemblies_fasta: {
          description: "Genomes to prepare for Genbank submission. One file per genome: all segments/chromosomes included in one file. All fasta files must contain exactly the same number of sequences as reference_fasta (which must equal the number of files in reference_annot_tbl).",
          patterns: ["*.fasta"]
        }
        reference_fastas: {
          description: "Reference genome, each segment/chromosome in a separate fasta file, in the exact same count and order as the segments/chromosomes described in genome_fasta. Headers must be Genbank accessions.",
          patterns: ["*.fasta"]
        }
        reference_feature_tables: {
          description: "NCBI Genbank feature table, each segment/chromosome in a separate TBL file, in the exact same count and order as the segments/chromosomes described in genome_fasta and reference_fastas. Accession numbers in the TBL files must correspond exactly to those in reference_fasta.",
          patterns: ["*.tbl"]
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
        taxid: {
          description: "The NCBI taxonomy ID for the species being submitted in this batch (all sequences in this batch must belong to the same taxid). https://www.ncbi.nlm.nih.gov/taxonomy/"
        }
        coverage_table: {
          description: "A two column tab text file mapping sample IDs (first column) to average sequencing coverage (second column, floating point number).",
          patterns: ["*.txt", "*.tsv"],
          category: "common"
        }
        sequencingTech: {
          description: "The type of sequencer used to generate reads. NCBI has a controlled vocabulary for this value which can be found here: https://submit.ncbi.nlm.nih.gov/structcomment/nongenomes/",
          category: "common"
        }
        organism: {
          description: "The scientific name for the organism being submitted. This is typically the species name and should match the name given by the NCBI Taxonomy database. For more info, see: https://www.ncbi.nlm.nih.gov/Sequin/sequin.hlp.html#Organism",
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

    call ncbi.biosample_to_genbank {
        input:
            biosample_attributes = biosample_attributes,
            num_segments         = length(reference_fastas),
            taxid                = taxid
    }

    scatter(assembly in assemblies_fasta) {
        call ncbi.align_and_annot_transfer_single as annot {
            input:
                genome_fasta             = assembly,
                reference_fastas         = reference_fastas,
                reference_feature_tables = reference_feature_tables
        }
    }

    call ncbi.generate_author_sbt_file as generate_author_sbt {
        input:
            author_list   = author_list,
            defaults_yaml = author_sbt_defaults_yaml,
            j2_template   = author_sbt_j2_template
    }

    call ncbi.prepare_genbank as prep_genbank {
        input:
            assemblies_fasta   = assemblies_fasta,
            annotations_tbl    = flatten(annot.genome_per_chr_tbls),
            authors_sbt        = generate_author_sbt.sbt_file,
            biosampleMap       = biosample_to_genbank.biosample_map,
            genbankSourceTable = biosample_to_genbank.genbank_source_modifier_table,
            coverage_table     = coverage_table,
            sequencingTech     = sequencingTech,
            comment            = comment,
            organism           = organism,
            molType            = molType
    }

    output {
        File        submission_zip         = prep_genbank.submission_zip
        File        archive_zip            = prep_genbank.archive_zip
        File        errorSummary           = prep_genbank.errorSummary
        
        File        biosample_map          = biosample_to_genbank.biosample_map
        File        genbank_source_table   = biosample_to_genbank.genbank_source_modifier_table
        
        Array[File] transferred_annot_tbls = flatten(annot.genome_per_chr_tbls)
        Array[File] genbank_preview_files  = prep_genbank.genbank_preview_files
        Array[File] validation_files       = prep_genbank.validation_files
        
        String      viral_phylo_version    = prep_genbank.viralngs_version
    }

}
