version 1.0

import "../tasks/tasks_ncbi.wdl" as ncbi
import "../tasks/tasks_nextstrain.wdl" as nextstrain

workflow genbank {

    meta {
        description: "Prepare assemblies for Genbank submission. This includes annotation by simple coordinate transfer from Genbank annotations and a multiple alignment. See https://viral-pipelines.readthedocs.io/en/latest/ncbi_submission.html for details."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
    }

    input {
        Array[File]+  assemblies_fasta

        File          authors_sbt
        File          biosample_attributes
        File          structured_comment_table
        Int           taxid=2697049
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
        authors_sbt: {
          description: "A genbank submission template file (SBT) with the author list, created at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/",
          patterns: ["*.sbt"]
        }
        biosample_attributes: {
          description: "A post-submission attributes file from NCBI BioSample, which is available at https://submit.ncbi.nlm.nih.gov/subs/ and clicking on 'Download attributes file with BioSample accessions'.",
          patterns: ["*.txt", "*.tsv"]
        }
        structured_comment_table: {
          description: "A six column tab text file with one row per sequence and the following header columns: SeqID, StructuredCommentPrefix, Assembly Method, Coverage        Sequencing Technology, StructuredCommentSuffix",
          patterns: ["*.txt", "*.tsv"],
          category: "common"
        }
        sequencingTech: {
          description: "The type of sequencer used to generate reads. NCBI has a controlled vocabulary for this value which can be found here: https://submit.ncbi.nlm.nih.gov/structcomment/nongenomes/",
          category: "common"
        }

    }

    scatter(assembly in assemblies_fasta) {
        call ncbi.vadr {
            input:
                genome_fasta = assembly
        }
        if (vadr.num_alerts==0) {
          File good_assembly = assembly
        }
    }

    call nextstrain.concatenate {
      input:
        infiles = select_all(good_assembly),
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

    call ncbi.package_genbank_ftp_submission {
      input:
        sequences_fasta = concatenate.combined,
        source_modifier_table = biosample_to_genbank.genbank_source_modifier_table,
        author_template_sbt = authors_sbt,
        structured_comment_table = structured_comment_table
    }

    output {
        File submission_zip = package_genbank_ftp_submission.submission_zip
        File submission_xml    = package_genbank_ftp_submission.submission_xml
        File submit_ready   = package_genbank_ftp_submission.submit_ready

        Int  num_successful = length(select_all(good_assembly))
        Int  num_input = length(assemblies_fasta)

        Array[File] vadr_outputs = vadr.outputs_tgz

        File biosample_map = biosample_to_genbank.biosample_map
        File genbank_source_table = biosample_to_genbank.genbank_source_modifier_table
    }

}
