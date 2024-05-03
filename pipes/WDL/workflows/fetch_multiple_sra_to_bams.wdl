version 1.0

import "../tasks/tasks_ncbi_tools.wdl" as ncbi_tools
import "../tasks/tasks_utils.wdl" as utils

workflow fetch_multiple_sra_to_bams {
    meta {
        description: "Retrieve reads for multiple SRA run IDs from the NCBI Short Read Archive in unaligned BAM format (multiple bam files) with relevant metadata encoded."
        author: "Broad Viral Genomics"
        email:  "viral-ngs@broadinstitute.org"
        allowNestedInputs: true
    }

    input {
        Array[String]+ SRA_IDs
    }

    parameter_meta {
        SRA_IDs: {
            description: "SRA run accessions (ex. *RR#######), NOT SRA study or sample accessions."
        }
    }

    scatter(sra_id in SRA_IDs) {
        call ncbi_tools.Fetch_SRA_to_BAM as scattered_fetch_sra_to_bam {
            input:SRA_ID = sra_id
        }

        Map[String,String] sra_outputs_map = {
            "reads_ubam":                scattered_fetch_sra_to_bam.reads_ubam,
            "sequencing_center":         scattered_fetch_sra_to_bam.sequencing_center,
            "sequencing_platform":       scattered_fetch_sra_to_bam.sequencing_platform,
            "sequencing_platform_model": scattered_fetch_sra_to_bam.sequencing_platform_model,
            "biosample_accession":       scattered_fetch_sra_to_bam.biosample_accession,
            "library_id":                scattered_fetch_sra_to_bam.library_id,
            "run_date":                  scattered_fetch_sra_to_bam.run_date,
            "sample_collection_date":    scattered_fetch_sra_to_bam.sample_collection_date,
            "sample_collected_by":       scattered_fetch_sra_to_bam.sample_collected_by,
            "sample_strain":             scattered_fetch_sra_to_bam.sample_strain,
            "sample_geo_loc":            scattered_fetch_sra_to_bam.sample_geo_loc,
            "sra_metadata":              scattered_fetch_sra_to_bam.sra_metadata
        }

        Array[String] metadata_for_accession = [
            sra_id,
            scattered_fetch_sra_to_bam.reads_ubam,
            scattered_fetch_sra_to_bam.sequencing_center,
            scattered_fetch_sra_to_bam.sequencing_platform,
            scattered_fetch_sra_to_bam.sequencing_platform_model,
            scattered_fetch_sra_to_bam.biosample_accession,
            scattered_fetch_sra_to_bam.library_id,
            scattered_fetch_sra_to_bam.run_date,
            scattered_fetch_sra_to_bam.sample_collection_date,
            scattered_fetch_sra_to_bam.sample_collected_by,
            scattered_fetch_sra_to_bam.sample_strain,
            scattered_fetch_sra_to_bam.sample_geo_loc,
            scattered_fetch_sra_to_bam.sra_metadata
        ]

        String sra_accession = sra_id
    }

    # create mapping from input SRA_ID to corresponding map of k:v containing metadata
    scatter(paired_metadata in zip(sra_accession, sra_outputs_map)){
        Map[String,Map[String,String]] combined_output_map = {
            paired_metadata.left: paired_metadata.right
        }
    }

    Array[String] metadata_header = [
                        "sra_run_accession",
                        "reads_ubam",
                        "sequencing_center",
                        "sequencing_platform",
                        "sequencing_platform_model",
                        "biosample_accession",
                        "library_id",
                        "run_date",
                        "sample_collection_date",
                        "sample_collected_by",
                        "sample_strain",
                        "sample_geo_loc",
                        "sra_metadata"
                    ]

    String input_ids_string = sep('_',SRA_IDs)

    call utils.concatenate as combined_metadata {
      input:
        # note that metadata_for_accession has type Array[Array[String]] since it is plural gathered scatter output
        infiles     = [write_tsv([metadata_header]), write_tsv(metadata_for_accession)],
        output_name = "run_metadata-${input_ids_string}.tsv"
    }

    output {
        # bam files for requested SRA IDs
        Array[File] read_bams = scattered_fetch_sra_to_bam.reads_ubam

        Array[ Map[ String, Map[String,String] ] ] collected_sra_metadata = combined_output_map
        File collected_sra_metadata_tsv = combined_metadata.combined
    }
}
