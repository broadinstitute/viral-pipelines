version 1.0

task Fetch_SRA_to_BAM {

    input {
        String  SRA_ID

        Int?    machine_mem_gb
        String  docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.1"
    }

    command {
        # pull reads from SRA and make a fully annotated BAM -- must succeed
        set -ex
        /opt/docker/scripts/sra_to_ubam.sh "~{SRA_ID}" "~{SRA_ID}.bam"

        # pull most metadata from BAM header
        set +e
        samtools view -H "~{SRA_ID}.bam" | grep ^@RG | head -1 | tr '\t' '\n' > header.txt
        grep CN header.txt | cut -f 2- -d : | tee OUT_CENTER
        grep PL header.txt | cut -f 2- -d : | tee OUT_PLATFORM
        grep SM header.txt | cut -f 2- -d : | tee OUT_BIOSAMPLE
        grep LB header.txt | cut -f 2- -d : | tee OUT_LIBRARY
        grep DT header.txt | cut -f 2 -d : | cut -f 1 -d T | tee OUT_RUNDATE
        samtools view -c "~{SRA_ID}.bam" | tee OUT_NUM_READS

        # pull other metadata from SRA -- allow for silent failures here!
        touch OUT_MODEL OUT_COLLECTION_DATE OUT_STRAIN OUT_COLLECTED_BY OUT_GEO_LOC
        esearch -db sra -q "~{SRA_ID}" | efetch -mode json -json > SRA.json
        cp SRA.json "~{SRA_ID}.json"
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
        jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.DESIGN.LIBRARY_DESCRIPTOR.LIBRARY_STRATEGY' \
            SRA.json | tee OUT_LIBRARY_STRATEGY

        python3 << CODE
        import json
        with open('SRA.json', 'rt') as inf:
            meta = json.load(inf)
            biosample = dict((x['TAG'],x['VALUE']) for x in meta['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['SAMPLE']['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE'])
            biosample['accession'] = meta['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['SAMPLE']['IDENTIFIERS']['EXTERNAL_ID']['content']
            biosample['message'] = 'Successfully loaded'
            biosample['bioproject_accession'] = meta['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['SAMPLE']['SAMPLE_LINKS']['SAMPLE_LINK']['XREF_LINK']['LABEL']
            biosample['sample_name'] = biosample['isolate']
            with open('~{SRA_ID}-biosample_attributes.json', 'wt') as outf:
                json.dump(biosample, outf)
        CODE
    }

    output {
        File    reads_ubam = "~{SRA_ID}.bam"
        Int     num_reads = read_int("OUT_NUM_READS")
        String  sequencing_center = read_string("OUT_CENTER")
        String  sequencing_platform = read_string("OUT_PLATFORM")
        String  sequencing_platform_model = read_string("OUT_MODEL")
        String  biosample_accession = read_string("OUT_BIOSAMPLE")
        String  library_id = read_string("OUT_LIBRARY")
        String  library_strategy = read_string("OUT_LIBRARY_STRATEGY")
        String  run_date = read_string("OUT_RUNDATE")
        String  sample_collection_date = read_string("OUT_COLLECTION_DATE")
        String  sample_collected_by = read_string("OUT_COLLECTED_BY")
        String  sample_strain = read_string("OUT_STRAIN")
        String  sample_geo_loc = read_string("OUT_GEO_LOC")
        File    sra_metadata = "~{SRA_ID}.json"
        File    biosample_attributes_json = "~{SRA_ID}-biosample_attributes.json"
    }

    runtime {
        cpu:     2
        memory:  select_first([machine_mem_gb, 6]) + " GB"
        disks:   "local-disk 750 LOCAL"
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
    }
}


task group_sra_bams_by_biosample {
  input {
    Array[File]   bam_filepaths
    Array[String] biosamples
    Array[File]   biosample_attributes_jsons
    Array[String] library_strategies
  }
  parameter_meta {
    bam_filepaths: {
      description: "all bam files",
      localization_optional: true,
      stream: true,
      patterns: ["*.bam"]
    }
  }
  command <<<
    python3 << CODE
    import os.path
    import json

    # WDL arrays to python arrays
    bam_uris = '~{sep="*" bam_filepaths}'.split('*')
    biosample_accs = '~{sep="*" biosamples}'.split('*')
    attributes = '~{sep="*" biosample_attributes_jsons}'.split('*')
    libstrats = '~{sep="*" library_strategies}'.split('*')
    assert len(bam_uris) == len(biosample_accs) == len(attributes) == len(libstrats)

    # lookup table files to dicts
    sample_to_bams = {}
    sample_to_attributes = {}
    sample_to_libstrat = {}
    attr_keys = set()
    for samn,bam,attr,libstrat in zip(biosample_accs,bam_uris, attributes, libstrats):
      sample_to_bams.setdefault(samn, [])
      sample_to_bams[samn].append(bam)
      sample_to_attributes[samn] = attr
      attr_keys.update(k for k,v in attr.items() if v)
      sample_to_libstrat.setdefault(samn, set())
      sample_to_libstrat[samn].add(libstrat)

    # write outputs
    with open('attributes.json') as out_attr:
        json.dump(sample_to_attributes, out_attr)
    with open('attributes.tsv') as out_attr:
        headers = tuple(sorted(attr_keys))
        out_attr.write('\t'.join(headers)+'\n')
        for sample in sorted(sample_to_bams.keys()):
            out_attr.write('\t'.join(sample_to_attributes[sample].get(h,'') for h in headers)+'\n')
    with open('grouped_bams', 'wt') as out_groups:
      with open('samns', 'wt') as out_samples:
        for sample in sorted(sample_to_bams.keys()):
          out_samples.write(sample+'\n')
          out_groups.write('\t'.join(sample_to_bams[sample])+'\n')
    with open('library_strategies.json') as out_attr:
        for k,v in sample_to_libstrat.items():
            sample_to_libstrat[k] = ';'.join(sorted(v))
        json.dump(sample_to_libstrat, out_attr)
    CODE
  >>>
  output {
    Array[Array[File]+] grouped_bam_filepaths = read_tsv('grouped_bams')
    Array[String]       biosample_accessions = read_lines('samns')
    Map[String,Map[String,String]] samn_to_attributes = read_json('attributes.json')
    File                biosample_attributes_tsv = 'attributes.tsv'
    Map[String,String]  samn_to_library_strategy = read_json('library_strategies.json')
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 100 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}

