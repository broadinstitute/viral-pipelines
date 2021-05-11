version 1.0

task Fetch_SRA_to_BAM {

    input {
        String  SRA_ID

        Int?    machine_mem_gb
        String  docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.6"
    }

    command <<<
        set -e
        # fetch SRA metadata on this record
        esearch -db sra -q "~{SRA_ID}" | efetch -mode json -json > SRA.json
        cp SRA.json "~{SRA_ID}.json"

        # pull reads from SRA and make a fully annotated BAM -- must succeed
        CENTER=$(jq -r .EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SUBMISSION.center_name SRA.json)
        PLATFORM=$(jq -r '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.PLATFORM | keys[] as $k | "\($k)"' SRA.json)
        MODEL=$(jq -r ".EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.PLATFORM.$PLATFORM.INSTRUMENT_MODEL" SRA.json)
        SAMPLE=$(jq -r '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.IDENTIFIERS.EXTERNAL_ID|select(.namespace == "BioSample")|.content' SRA.json)
        LIBRARY=$(jq -r .EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.alias SRA.json)
        RUNDATE=$(jq -r '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.RUN_SET.RUN.SRAFiles|if (.SRAFile|type) == "object" then .SRAFile.date else [.SRAFile[]|select(.supertype == "Original")][0].date end' SRA.json | cut -f 1 -d ' ')

        if [ "$PLATFORM" = "OXFORD_NANOPORE" ]; then
            # per the SAM/BAM specification
            SAM_PLATFORM="ONT"
        else
            SAM_PLATFORM="$PLATFORM"
        fi

        sam-dump --unaligned --header "~{SRA_ID}" \
            | samtools view -bhS - \
            > temp.bam
        picard AddOrReplaceReadGroups \
            I=temp.bam \
            O="~{SRA_ID}.bam" \
            RGID=1 \
            RGLB="$LIBRARY" \
            RGSM="$SAMPLE" \
            RGPL="$SAM_PLATFORM" \
            RGPU="$LIBRARY" \
            RGPM="$MODEL" \
            RGDT="$RUNDATE" \
            RGCN="$CENTER" \
            VALIDATION_STRINGENCY=SILENT
        rm temp.bam
        samtools view -H "~{SRA_ID}.bam"

        # emit numeric WDL outputs
        echo $CENTER > OUT_CENTER
        echo $PLATFORM > OUT_PLATFORM
        echo $SAMPLE > OUT_BIOSAMPLE
        echo $LIBRARY > OUT_LIBRARY
        echo $RUNDATE > OUT_RUNDATE
        samtools view -c "~{SRA_ID}.bam" | tee OUT_NUM_READS

        # pull other metadata from SRA -- allow for silent failures here!
        touch OUT_MODEL OUT_COLLECTION_DATE OUT_STRAIN OUT_COLLECTED_BY OUT_GEO_LOC
        set +e
        jq -r \
            .EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.PLATFORM."$PLATFORM".INSTRUMENT_MODEL \
            SRA.json | tee OUT_MODEL
        jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "collection_date" or .TAG=="collection date")|.VALUE' \
            SRA.json | tee OUT_COLLECTION_DATE
        jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "strain")|.VALUE' \
            SRA.json | tee OUT_STRAIN
        jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "collected_by" or .TAG == "collecting institution")|.VALUE' \
            SRA.json | tee OUT_COLLECTED_BY
        jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE[]|select(.TAG == "geo_loc_name" or .TAG == "geographic location (country and/or sea)")|.VALUE' \
            SRA.json | tee OUT_GEO_LOC
        jq -r \
            '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.DESIGN.LIBRARY_DESCRIPTOR.LIBRARY_STRATEGY' \
            SRA.json | tee OUT_LIBRARY_STRATEGY

        set -e
        python3 << CODE
        import json
        with open('SRA.json', 'rt') as inf:
            meta = json.load(inf)
        # reorganize to look more like a biosample attributes tsv
        biosample = dict((x['TAG'],x['VALUE']) for x in meta['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['SAMPLE']['SAMPLE_ATTRIBUTES']['SAMPLE_ATTRIBUTE'])
        biosample['accession'] = meta['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['SAMPLE']['IDENTIFIERS']['EXTERNAL_ID']['content']
        biosample['message'] = 'Successfully loaded'
        biosample['bioproject_accession'] = meta['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['STUDY']['IDENTIFIERS']['EXTERNAL_ID']['content']
        biosample['sample_name'] = biosample['isolate']
        for k,v in biosample.items():
            if v == 'not provided':
                biosample[k] = ''

        # British to American conversions (NCBI vs ENA)
        us_to_uk = {
            'sample_name': 'Sample Name',
            'isolate': 'Sample Name',
            'collected_by': 'collecting institution',
            'collection_date': 'collection date',
            'geo_loc_name': 'geographic location (country and/or sea)',
            'host': 'host scientific name',
        }
        for key_us, key_uk in us_to_uk.items():
            if not biosample.get(key_us,''):
                biosample[key_us] = biosample.get(key_uk,'')

        # write outputs
        with open('~{SRA_ID}-biosample_attributes.json', 'wt') as outf:
            json.dump(biosample, outf)
        CODE
    >>>

    output {
        File    reads_ubam                = "~{SRA_ID}.bam"
        Int     num_reads                 = read_int("OUT_NUM_READS")
        String  sequencing_center         = read_string("OUT_CENTER")
        String  sequencing_platform       = read_string("OUT_PLATFORM")
        String  sequencing_platform_model = read_string("OUT_MODEL")
        String  biosample_accession       = read_string("OUT_BIOSAMPLE")
        String  library_id                = read_string("OUT_LIBRARY")
        String  library_strategy          = read_string("OUT_LIBRARY_STRATEGY")
        String  run_date                  = read_string("OUT_RUNDATE")
        String  sample_collection_date    = read_string("OUT_COLLECTION_DATE")
        String  sample_collected_by       = read_string("OUT_COLLECTED_BY")
        String  sample_strain             = read_string("OUT_STRAIN")
        String  sample_geo_loc            = read_string("OUT_GEO_LOC")
        File    sra_metadata              = "~{SRA_ID}.json"
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

task fetch_biosamples {

    input {
        Array[String]  biosample_ids

        String         out_basename = "biosample_attributes"
        String         docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.6"
    }

    command <<<
        set -e
        /opt/docker/scripts/biosample-fetch_attributes.py \
            ~{sep=' ' biosample_ids} "~{out_basename}"
    >>>

    output {
        File    biosample_attributes_tsv  = "~{out_basename}.tsv"
        File    biosample_attributes_json = "~{out_basename}.json"
    }

    runtime {
        cpu:     2
        memory:  "3 GB"
        disks:   "local-disk 50 HDD"
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
    }
}

task ncbi_ftp_upload {
    input {
        Array[File]    submit_files
        File           config_js
        String         target_path

        String         wait_for="1"  # all, disabled, some number

        String         docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.6"
    }

    command <<<
        set -e
        cd /opt/converter
        cp "~{config_js}" src/config.js
        rm -f files/sample.tsv reports/sample-report.xml
        cp ~{sep=' ' submit_files} files/
        MANIFEST=$(ls -1 files | paste -sd,)
        echo "uploading: $MANIFEST to destination ftp folder ~{target_path}"
        echo "Asymmetrik script version: $ASYMMETRIK_REPO_COMMIT"
        node src/main.js --debug \
            --uploadFiles="$MANIFEST" \
            --poll="~{wait_for}" \
            --uploadFolder="~{target_path}"
        cd -
        cp /opt/converter/reports/*report*.xml .
        ls -alF files reports
    >>>

    output {
        Array[File] reports_xmls = glob("*report*.xml")
    }

    runtime {
        cpu:     2
        memory:  "2 GB"
        disks:   "local-disk 100 HDD"
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
    }
}

task biosample_submit_tsv_to_xml {
    input {
        File     meta_submit_tsv
        File     config_js

        String   docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.6"
    }
    command <<<
        set -e
        cd /opt/converter
        cp "~{config_js}" src/config.js
        rm files/sample.tsv
        cp "~{meta_submit_tsv}" files/
        echo "Asymmetrik script version: $ASYMMETRIK_REPO_COMMIT"
        node src/main.js --debug \
            -i=$(basename "~{meta_submit_tsv}") \
            --runTestMode=true
        cd -
        cp "/opt/converter/files/~{basename(meta_submit_tsv, '.tsv')}-submission.xml" .
    >>>
    output {
        File   submission_xml = "~{basename(meta_submit_tsv, '.tsv')}-submission.xml"
    }
    runtime {
        cpu:     2
        memory:  "2 GB"
        disks:   "local-disk 100 HDD"
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
    }
}

task biosample_submit_tsv_ftp_upload {
    input {
        File     meta_submit_tsv
        File     config_js
        String   target_path

        String   docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.6"
    }
    String base=basename(meta_submit_tsv, '.tsv')
    command <<<
        set -e
        cd /opt/converter
        cp "~{config_js}" src/config.js
        rm -f files/sample.tsv reports/sample-report.xml
        cp "~{meta_submit_tsv}" files/
        echo "Asymmetrik script version: $ASYMMETRIK_REPO_COMMIT"
        node src/main.js --debug \
            -i=$(basename "~{meta_submit_tsv}") \
            --uploadFolder="~{target_path}"
        cd -
        cp /opt/converter/reports/~{base}-attributes.tsv /opt/converter/files/~{base}-submission.xml /opt/converter/reports/~{base}-report.*.xml .
    >>>
    output {
        File        attributes_tsv = "~{base}-attributes.tsv"
        File        submission_xml = "~{base}-submission.xml"
        Array[File] reports_xmls   = glob("~{base}-report.*.xml")
    }
    runtime {
        cpu:     2
        memory:  "2 GB"
        disks:   "local-disk 100 HDD"
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
    }
}

task biosample_xml_response_to_tsv {
    input {
        File     meta_submit_tsv
        File     ncbi_report_xml

        String   docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.6"
    }
    String out_name = "~{basename(meta_submit_tsv, '.tsv')}-attributes.tsv"
    command <<<
        set -e
        cd /opt/converter
        cp "~{meta_submit_tsv}" files/submit.tsv
        cp "~{ncbi_report_xml}" reports/report.xml
        echo "Asymmetrik script version: $ASYMMETRIK_REPO_COMMIT"
        node src/main.js --debug \
            -i=submit.tsv \
            -p=report.xml
        cd -
        cp /opt/converter/reports/submit-attributes.tsv "~{out_name}"
    >>>
    output {
        File   biosample_attributes_tsv = "~{out_name}"
    }
    runtime {
        cpu:     2
        memory:  "2 GB"
        disks:   "local-disk 100 HDD"
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
    Array[String] seq_platforms
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
    bam_uris = list(x for x in '~{sep="*" bam_filepaths}'.split('*') if x)
    biosample_accs = list(x for x in '~{sep="*" biosamples}'.split('*') if x)
    attributes = list(x for x in '~{sep="*" biosample_attributes_jsons}'.split('*') if x)
    libstrats = list(x for x in '~{sep="*" library_strategies}'.split('*') if x)
    seqplats = list(x for x in '~{sep="*" seq_platforms}'.split('*') if x)
    assert len(bam_uris) == len(biosample_accs) == len(attributes) == len(libstrats) == len(seqplats)

    # lookup table files to dicts
    sample_to_bams = {}
    sample_to_attributes = {}
    sample_to_libstrat = {}
    sample_to_seqplat = {}
    attr_keys = set()
    for samn,bam,attr_file,libstrat,seqplat in zip(biosample_accs,bam_uris, attributes, libstrats, seqplats):
      sample_to_bams.setdefault(samn, [])
      sample_to_bams[samn].append(bam)
      with open(attr_file, 'rt') as inf:
        attr = json.load(inf)
      sample_to_attributes[samn] = attr
      attr_keys.update(k for k,v in attr.items() if v)
      sample_to_libstrat.setdefault(samn, set())
      sample_to_libstrat[samn].add(libstrat)
      sample_to_seqplat.setdefault(samn, set())
      sample_to_seqplat[samn].add(seqplat)

    # write outputs
    with open('attributes.json', 'wt') as outf:
        json.dump(sample_to_attributes, outf)
    with open('attributes.tsv', 'wt') as outf:
        headers = tuple(sorted(attr_keys))
        outf.write('\t'.join(headers)+'\n')
        for sample in sorted(sample_to_bams.keys()):
            outf.write('\t'.join(sample_to_attributes[sample].get(h,'') for h in headers)+'\n')
    with open('grouped_bams', 'wt') as out_groups:
      with open('samns', 'wt') as out_samples:
        for sample in sorted(sample_to_bams.keys()):
          out_samples.write(sample+'\n')
          out_groups.write('\t'.join(sample_to_bams[sample])+'\n')
    with open('library_strategies.json', 'wt') as outf:
        for k,v in sample_to_libstrat.items():
            sample_to_libstrat[k] = ';'.join(sorted(v))
        json.dump(sample_to_libstrat, outf)
    with open('sequencing_platforms.json', 'wt') as outf:
        for k,v in sample_to_seqplat.items():
            sample_to_seqplat[k] = ';'.join(sorted(v))
        json.dump(sample_to_seqplat, outf)
    CODE
  >>>
  output {
    Array[Array[File]+]            grouped_bam_filepaths       = read_tsv('grouped_bams')
    Array[String]                  biosample_accessions        = read_lines('samns')
    Map[String,Map[String,String]] samn_to_attributes          = read_json('attributes.json')
    File                           biosample_attributes_tsv    = 'attributes.tsv'
    Map[String,String]             samn_to_library_strategy    = read_json('library_strategies.json')
    Map[String,String]             samn_to_sequencing_platform = read_json('sequencing_platforms.json')
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 100 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}

