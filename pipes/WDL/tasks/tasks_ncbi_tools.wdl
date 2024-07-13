version 1.0

task Fetch_SRA_to_BAM {

    input {
        String  SRA_ID

        String? sample_name
        String? email_address
        String? ncbi_api_key
        Int?    machine_mem_gb
        String  docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.10"
    }
    Int disk_size = 750
    meta {
        description: "This searches NCBI SRA for accessions using the Entrez interface, collects associated metadata, and returns read sets as unaligned BAM files with metadata loaded in. Useful metadata from BioSample is also output from this task directly. This has been tested with both SRA and ENA accessions. This queries the NCBI production database, and as such, the output of this task is non-deterministic given the same input."
        volatile: true
    }
    command <<<
        set -e
        ~{if defined(ncbi_api_key) then "export NCBI_API_KEY=~{ncbi_api_key}}" else ""}

        # fetch SRA metadata on this record
        esearch ~{if defined(email_address) then "-email ~{email_address}" else ""} -db sra -query "~{SRA_ID}" | efetch -db sra ~{if defined(email_address) then "-email ~{email_address}" else ""} -mode json -json > SRA.json

        cp SRA.json "~{SRA_ID}.json"

        # pull reads from SRA and make a fully annotated BAM -- must succeed
        CENTER=$(jq -r .EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SUBMISSION.center_name SRA.json)
        PLATFORM=$(jq -r '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.PLATFORM | keys[] as $k | "\($k)"' SRA.json)
        MODEL=$(jq -r ".EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.PLATFORM.$PLATFORM.INSTRUMENT_MODEL" SRA.json)
        SAMPLE=$(jq -r '.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.SAMPLE.IDENTIFIERS.EXTERNAL_ID[]|select(.namespace == "BioSample")|.content' SRA.json)
        LIBRARY=$(jq -r .EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.EXPERIMENT.alias SRA.json)
        RUNDATE=$(jq -r '(.EXPERIMENT_PACKAGE_SET.EXPERIMENT_PACKAGE.RUN_SET | (if (.RUN|type) == "object" then (.RUN) else (.RUN[] | select(any(.; .accession == "~{SRA_ID}"))) end) | .SRAFiles) | if (.SRAFile|type) == "object" then .SRAFile.date else [.SRAFile[]|select(.supertype == "Original" or .supertype=="Primary ETL")][0].date end' SRA.json | cut -f 1 -d ' ')
        
    >>>

    output {
        File    reads_ubam                = "~{SRA_ID}.bam"
    }

    runtime {
        cpu:     2
        memory:  select_first([machine_mem_gb, 6]) + " GB"
        disks:   "local-disk " + disk_size + " LOCAL"
        disk:    disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
        maxRetries: 2
    }
}

task fetch_genbank_metadata {
    input {
        String  genbank_accession
        String  docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.10"
    }
    Int disk_size = 50
    command <<<
        set -e
        source /opt/miniconda/bin/activate $CONDA_DEFAULT_ENV # for miniwdl / non-login docker runners
        esearch -db nuccore -q "~{genbank_accession}" | efetch -db nuccore -format gb -mode xml -json  > gb.json
        jq -r '[.GBSet.GBSeq."GBSeq_feature-table".GBFeature[0].GBFeature_quals.GBQualifier|.[]|{(.GBQualifier_name): .GBQualifier_value}]|add ' gb.json > "~{genbank_accession}".metadata.json
        jq -r '.db_xref' "~{genbank_accession}".metadata.json | grep ^taxon: | cut -f 2 -d : > taxid.txt
        jq -r '.organism' "~{genbank_accession}".metadata.json > organism.txt
    >>>
    output {
        Map[String,String] metadata = read_json("~{genbank_accession}.metadata.json")
        String taxid = read_string("taxid.txt")
        String organism = read_string("organism.txt")
    }
    runtime {
        cpu:     1
        memory:  "1 GB"
        disks:   "local-disk " + disk_size + " LOCAL"
        disk:    disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        docker:  docker
        maxRetries: 2
    }
}

task biosample_tsv_filter_preexisting {
    input {
        File           meta_submit_tsv

        String         out_basename = "biosample_attributes"
        String         docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.10"
    }
    Int disk_size = 50
    meta {
        description: "This task takes a metadata TSV for submission to NCBI BioSample and filters out entries that have already been submitted to NCBI. This queries the NCBI production database, and as such, the output of this task is non-deterministic given the same input."
        volatile: true
    }    
    command <<<
        set -e

        # extract all the sample_name values from input tsv
        python3<<CODE
        import csv
        with open('~{meta_submit_tsv}', 'rt') as inf:
            with open('SAMPLES.txt', 'w', newline='') as outf:
                for row in csv.DictReader(inf, delimiter='\t'):
                    outf.write(row['isolate'] + '\n')
        CODE
        cat SAMPLES.txt | wc -l | tee COUNT_IN

        # fetch attributes file for anything already registered
        /opt/docker/scripts/biosample-fetch_attributes.py \
            $(cat SAMPLES.txt) "~{out_basename}"

        # extract all the sample_name values from NCBI
        python3<<CODE
        import csv
        with open('~{out_basename}.tsv', 'rt') as inf:
            with open('FOUND.txt', 'w', newline='') as outf:
                for row in csv.DictReader(inf, delimiter='\t'):
                    outf.write(row['isolate'] + '\n')
        CODE
        cat FOUND.txt | wc -l | tee COUNT_FOUND

        # filter out from input
        set +e
        grep -v -F -f FOUND.txt "~{meta_submit_tsv}" > "~{out_basename}.unsubmitted.tsv"
        cat "~{out_basename}.unsubmitted.tsv" | wc -l | tee COUNT_UNFOUND
    >>>
    output {
        File    meta_unsubmitted_tsv = "~{out_basename}.unsubmitted.tsv"
        File    biosample_attributes_tsv  = "~{out_basename}.tsv"
        Int     num_in = read_int("COUNT_IN")
        Int     num_found = read_int("COUNT_FOUND")
        Int     num_not_found = read_int("COUNT_UNFOUND") - 1
    }
    runtime {
        cpu:     2
        memory:  "3 GB"
        disks:   "local-disk " + disk_size + " HDD"
        disk:    disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
        maxRetries: 2
    }
}

task fetch_biosamples {
    input {
        Array[String]  biosample_ids

        String         out_basename = "biosample_attributes"
        String         docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.10"
    }
    Int disk_size = 50
    meta {
        description: "This searches NCBI BioSample for accessions or keywords using the Entrez interface and returns any hits in the form of a BioSample attributes TSV. This queries the NCBI production database, and as such, the output of this task is non-deterministic given the same input."
        volatile: true
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
        disks:   "local-disk " + disk_size + " HDD"
        disk:   disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
        maxRetries: 2
    }
}

task ncbi_sftp_upload {
    input {
        File           submission_xml
        Array[File]    additional_files = []
        File           config_js
        String         target_path

        String         wait_for="1"  # all, disabled, some number

        String         docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.10"
    }
    Int disk_size = 100
    command <<<
        set -e
        cd /opt/converter
        cp "~{config_js}" src/config.js
        rm -rf files/tests
        cp "~{submission_xml}" files/submission.xml
        if [[ "~{length(additional_files)}" != "0" ]]; then
            cp ~{sep=' ' additional_files} files/
        fi
        MANIFEST=$(ls -1 files | paste -sd,)
        echo "uploading: $MANIFEST to destination ftp folder ~{target_path}"
        echo "Asymmetrik script version: $ASYMMETRIK_REPO_COMMIT"
        node src/main.js --debug \
            --uploadFiles="$MANIFEST" \
            --poll="~{wait_for}" \
            --uploadFolder="~{target_path}"
        ls -alF files reports
        cd -
        cp /opt/converter/reports/*report*.xml .
    >>>

    output {
        Array[File] reports_xmls = glob("*report*.xml")
    }

    runtime {
        cpu:     2
        memory:  "2 GB"
        disks:   "local-disk " + disk_size + " HDD"
        disk:    disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
        maxRetries: 0
    }
}

task sra_tsv_to_xml {
    input {
        File     meta_submit_tsv
        File     config_js
        String   bioproject
        String   data_bucket_uri

        String   docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.10"
    }
    Int disk_size = 50
    command <<<
        set -e
        cd /opt/converter
        cp "~{config_js}" src/config.js
        cp "~{meta_submit_tsv}" files/
        echo "Asymmetrik script version: $ASYMMETRIK_REPO_COMMIT"
        node src/main.js --debug \
            -i=$(basename "~{meta_submit_tsv}") \
            --submissionType=sra \
            --bioproject="~{bioproject}" \
            --submissionFileLoc="~{data_bucket_uri}" \
            --runTestMode=true
        cd -
        cp "/opt/converter/files/~{basename(meta_submit_tsv, '.tsv')}-submission.xml" .
    >>>
    output {
        File   submission_xml = "~{basename(meta_submit_tsv, '.tsv')}-submission.xml"
    }
    runtime {
        cpu:     1
        memory:  "2 GB"
        disks:   "local-disk " + disk_size + " HDD"
        disk:    disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
        maxRetries: 2
    }
}

task biosample_submit_tsv_to_xml {
    input {
        File     meta_submit_tsv
        File     config_js

        String   docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.10"
    }
    Int disk_size = 50
    meta {
        description: "This converts a web portal submission TSV for NCBI BioSample into an ftp-appropriate XML submission for NCBI BioSample. It does not connect to NCBI, and does not submit or fetch any data."
    }
    command <<<
        set -e
        cd /opt/converter
        cp "~{config_js}" src/config.js
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
        cpu:     1
        memory:  "2 GB"
        disks:   "local-disk " + disk_size + " HDD"
        disk:    disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
        maxRetries: 2
    }
}

task biosample_submit_tsv_ftp_upload {
    input {
        File     meta_submit_tsv
        File     config_js
        String   target_path

        String   docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.10"
    }
    String base=basename(meta_submit_tsv, '.tsv')
    Int disk_size = 100    
    meta {
        description: "This registers a table of metadata with NCBI BioSample. It accepts a TSV similar to the web UI input at submit.ncbi.nlm.nih.gov, but converts to an XML, submits via their FTP/XML API, awaits a response, and retrieves a resulting attributes table and returns that as a TSV. This task registers live data with the production NCBI database."
    }
    command <<<
        set -e
        cd /opt/converter
        cp "~{config_js}" src/config.js
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
        disks:   "local-disk " + disk_size + " HDD"
        disk:    disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
        maxRetries: 0
    }
}

task biosample_xml_response_to_tsv {
    input {
        File     meta_submit_tsv
        File     ncbi_report_xml

        String   docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.10"
    }
    String out_name = "~{basename(meta_submit_tsv, '.tsv')}-attributes.tsv"
    Int disk_size = 100
    meta {
        description: "This converts an FTP-based XML response from BioSample into a web-portal-style attributes.tsv file with metadata and accessions. This task does not communicate with NCBI, it only parses pre-retrieved responses."
    }
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
        disks:   "local-disk " + disk_size + " HDD"
        disk:    disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        docker:  docker
        maxRetries: 2
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
  Int disk_size = 100
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
    disks:   "local-disk " + disk_size + " HDD"
    disk:    disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

