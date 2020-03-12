version 1.0

task Fetch_SRA_to_BAM {

    input {
        String  SRA_ID
        String  docker = "quay.io/broadinstitute/ncbi-tools"
    }

    command {
        set -ex -o pipefail

        # pull reads from SRA and make a fully annotated BAM
        /opt/docker/scripts/sra_to_ubam.sh ${SRA_ID} ${SRA_ID}.bam

        # pull other metadata from SRA
        esearch -db sra -q "${SRA_ID}" | efetch -mode json -json > ${SRA_ID}.json

        cat ${SRA_ID}.json | jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SUBMISSION.center_name' \
            | tee OUT_CENTER
        cat ${SRA_ID}.json | jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.PLATFORM | keys[] as $k | "\($k)"' \
            | tee OUT_PLATFORM
        cat ${SRA_ID}.json | jq -r \
            .EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.PLATFORM.$PLATFORM.INSTRUMENT_MODEL \
            | tee OUT_MODEL
        cat ${SRA_ID}.json | jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.IDENTIFIERS.EXTERNAL_ID|select(.namespace == "BioSample")|.content' \
            | tee OUT_BIOSAMPLE
        cat ${SRA_ID}.json | jq -r \
            .EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.DESIGN.LIBRARY_DESCRIPTOR.LIBRARY_NAME \
            | tee OUT_LIBRARY
        cat ${SRA_ID}.json | jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.RUN_SET.RUN.SRAFiles.SRAFile[]|select(.supertype == "Original")|.date' \
            | cut -f 1 -d ' ' \
            | tee OUT_RUNDATE
        cat ${SRA_ID}.json | jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "collection_date")|.VALUE' \
            | tee OUT_COLLECTION_DATE
        cat ${SRA_ID}.json | jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "strain")|.VALUE' \
            | tee OUT_STRAIN
        cat ${SRA_ID}.json | jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "collected_by")|.VALUE' \
            | tee OUT_COLLECTED_BY
        cat ${SRA_ID}.json | jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "geo_loc_name")|.VALUE' \
            | tee OUT_GEO_LOC
    }

    output {
        File    reads_ubam = "${SRA_ID}.bam"
        String  sequencing_center = read_string("OUT_CENTER")
        String  sequencing_platform = read_string("OUT_PLATFORM")
        String  sequencing_platform_model = read_string("OUT_MODEL")
        String  biosample_accession = read_string("OUT_BIOSAMPLE")
        String  library_id = read_string("OUT_LIBRARY")
        String  run_date = read_string("OUT_RUNDATE")
        String  sample_collection_date = read_string("OUT_COLLECTION_DATE")
        String  sample_collected_by = read_string("OUT_COLLECTED_BY")
        String  sample_strain = read_string("OUT_STRAIN")
        String  sample_geo_loc = read_string("OUT_GEO_LOC")
        File    sra_metadata = "${SRA_ID}.json"
    }

    runtime {
        cpu:     4
        memory:  "15 GB"
        disks:   "local-disk 750 LOCAL"
        docker:  "${docker}"
    }
}
