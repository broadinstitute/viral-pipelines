version 1.0

task Fetch_SRA_to_BAM {

    input {
        String  SRA_ID

        Int?    machine_mem_gb
        String  docker = "quay.io/broadinstitute/ncbi-tools"
    }

    command {
        # pull reads from SRA and make a fully annotated BAM -- must succeed
        set -ex
        /opt/docker/scripts/sra_to_ubam.sh "${SRA_ID}" "${SRA_ID}.bam"

        # pull most metadata from BAM header
        set +e
        samtools view -H "${SRA_ID}.bam" | grep ^@RG | head -1 | tr '\t' '\n' > header.txt
        grep CN header.txt | cut -f 2- -d : | tee OUT_CENTER
        grep PL header.txt | cut -f 2- -d : | tee OUT_PLATFORM
        grep SM header.txt | cut -f 2- -d : | tee OUT_BIOSAMPLE
        grep LB header.txt | cut -f 2- -d : | tee OUT_LIBRARY
        grep DT header.txt | cut -f 2 -d : | cut -f 1 -d T | tee OUT_RUNDATE

        # pull other metadata from SRA -- allow for silent failures here!
        touch OUT_MODEL OUT_COLLECTION_DATE OUT_STRAIN OUT_COLLECTED_BY OUT_GEO_LOC
        esearch -db sra -q "${SRA_ID}" | efetch -mode json -json > SRA.json
        jq -r \
            .EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.PLATFORM."$(<OUT_PLATFORM)".INSTRUMENT_MODEL \
            SRA.json | tee OUT_MODEL
        jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "collection_date")|.VALUE' \
            SRA.json | tee OUT_COLLECTION_DATE
        jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "strain")|.VALUE' \
            SRA.json | tee OUT_STRAIN
        jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "collected_by")|.VALUE' \
            SRA.json | tee OUT_COLLECTED_BY
        jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "geo_loc_name")|.VALUE' \
            SRA.json | tee OUT_GEO_LOC
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
        cpu:     2
        memory:  select_first([machine_mem_gb, 6]) + " GB"
        disks:   "local-disk 750 LOCAL"
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  "${docker}"
    }
}
