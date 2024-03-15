version 1.0

task nextclade_one_sample {
    meta {
        description: "Nextclade classification of one sample. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
    }
    input {
        File  genome_fasta
        File? root_sequence
        File? auspice_reference_tree_json
        File? qc_config_json
        File? gene_annotations_json
        File? pcr_primers_csv
        File? virus_properties
        String? dataset_name
        Int    disk_size = 50
        String docker = "nextstrain/nextclade:2.14.0"
    }
    String basename = basename(genome_fasta, ".fasta")
    command {
        set -e
        apt-get update
        apt-get -y install python3

        export NEXTCLADE_VERSION="$(nextclade --version)"
        echo $NEXTCLADE_VERSION > VERSION

        # grab named nextclade dataset
        DATASET_ARG=""
        if [ -n "~{default='' dataset_name}" ]; then
            nextclade dataset get --name="~{default='' dataset_name}" --output-dir=.
            DATASET_ARG="--input-dataset ."
        python3<<CODE1
        import json, os
        with open('tag.json', 'rt') as inf:
            datasetinfo = json.load(inf)
        with open('VERSION', 'wt') as outf:
            outf.write(os.environ['NEXTCLADE_VERSION'] + "; name=" + datasetinfo['name'] + "; tag=" + datasetinfo['tag'] + "\n")
        CODE1
        fi

        nextclade run \
            $DATASET_ARG \
            ~{"--input-root-seq " + root_sequence} \
            ~{"--input-tree " + auspice_reference_tree_json} \
            ~{"--input-qc-config " + qc_config_json} \
            ~{"--input-gene-map " + gene_annotations_json} \
            ~{"--input-pcr-primers " + pcr_primers_csv} \
            ~{"--input-virus-properties " + virus_properties}  \
            --output-all=. \
            --output-basename "~{basename}" \
            --output-json "~{basename}".nextclade.json \
            --output-tsv  "~{basename}".nextclade.tsv \
            --output-tree "~{basename}".nextclade.auspice.json \
            "~{genome_fasta}"
        python3 <<CODE
        # transpose table
        import codecs
        with codecs.open('~{basename}.nextclade.tsv', 'r', encoding='utf-8') as inf:
            with codecs.open('transposed.tsv', 'w', encoding='utf-8') as outf:
                for c in zip(*(l.rstrip().split('\t') for l in inf)):
                    outf.write('\t'.join(c)+'\n')
        CODE
        grep ^clade\\W transposed.tsv | cut -f 2 | grep -v clade > NEXTCLADE_CLADE
        grep ^aaSubstitutions\\W transposed.tsv | cut -f 2 | grep -v aaSubstitutions > NEXTCLADE_AASUBS
        grep ^aaDeletions\\W transposed.tsv | cut -f 2 | grep -v aaDeletions > NEXTCLADE_AADELS
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu:    2
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        String nextclade_version = read_string("VERSION")
        File   nextclade_json    = "~{basename}.nextclade.json"
        File   auspice_json      = "~{basename}.nextclade.auspice.json"
        File   nextclade_tsv     = "~{basename}.nextclade.tsv"
        String nextclade_clade   = read_string("NEXTCLADE_CLADE")
        String aa_subs_csv       = read_string("NEXTCLADE_AASUBS")
        String aa_dels_csv       = read_string("NEXTCLADE_AADELS")
    }
}

task nextclade_many_samples {
    meta {
        description: "Nextclade classification of many samples. Leaving optional inputs unspecified will use SARS-CoV-2 defaults."
    }
    input {
        Array[File]+ genome_fastas
        File?        root_sequence
        File?        auspice_reference_tree_json
        File?        qc_config_json
        File?        gene_annotations_json
        File?        pcr_primers_csv
        File?        virus_properties
        String?      dataset_name
        String       basename
        File?        genome_ids_setdefault_blank
        Int          disk_size = 150
        String       docker = "nextstrain/nextclade:2.14.0"
    }
    command <<<
        set -e
        apt-get update
        apt-get -y install python3

        export NEXTCLADE_VERSION="$(nextclade --version)"
        echo $NEXTCLADE_VERSION > VERSION

        # grab named nextclade dataset
        DATASET_ARG=""
        if [ -n "~{default='' dataset_name}" ]; then
            nextclade dataset get --name="~{default='' dataset_name}" --output-dir=.
            DATASET_ARG="--input-dataset ."
        python3<<CODE1
        import json, os
        with open('tag.json', 'rt') as inf:
            datasetinfo = json.load(inf)
        with open('VERSION', 'wt') as outf:
            outf.write(os.environ['NEXTCLADE_VERSION'] + "; name=" + datasetinfo['name'] + "; tag=" + datasetinfo['tag'] + "\n")
        CODE1
        fi

        nextclade run \
            $DATASET_ARG \
            ~{"--input-root-seq " + root_sequence} \
            ~{"--input-tree " + auspice_reference_tree_json} \
            ~{"--input-qc-config " + qc_config_json} \
            ~{"--input-gene-map " + gene_annotations_json} \
            ~{"--input-pcr-primers " + pcr_primers_csv} \
            ~{"--input-virus-properties " + virus_properties}  \
            --output-all=. \
            --output-basename "~{basename}" \
            --output-json "~{basename}".nextclade.json \
            --output-tsv  "~{basename}".nextclade.tsv \
            --output-tree "~{basename}".nextclade.auspice.json \
            --output-fasta "~{basename}".nextalign.msa.fasta \
            ~{sep=" " genome_fastas}

        python3 <<CODE
        # transpose table
        import codecs, csv, json
        out_maps = {'clade':{}, 'aaSubstitutions':{}, 'aaDeletions':{}}
        with codecs.open('IDLIST', 'w', encoding='utf-8') as outf_ids:
            # parse entries from output tsv into jsons
            with codecs.open('~{basename}.nextclade.tsv', 'r', encoding='utf-8') as inf:
                for row in csv.DictReader(inf, delimiter='\t'):
                    for k in out_maps.keys():
                        out_maps[k][row['seqName']] = row[k]
                    outf_ids.write(row['seqName']+'\n')
            # fill empty values for anything not mentioned by output tsv
            with codecs.open("~{default='/dev/null' genome_ids_setdefault_blank}", 'r', encoding='utf-8') as inf:
                for line in inf:
                    seqName = line.strip()
                    if seqName and (seqName not in out_maps['clade']):
                        for k in out_maps.keys():
                            out_maps[k][seqName] = ''
                        outf_ids.write(seqName+'\n')

        with codecs.open('NEXTCLADE_CLADE.json', 'w', encoding='utf-8') as outf:
            json.dump(out_maps['clade'], outf)
        with codecs.open('NEXTCLADE_AASUBS.json', 'w', encoding='utf-8') as outf:
            json.dump(out_maps['aaSubstitutions'], outf)
        with codecs.open('NEXTCLADE_AADELS.json', 'w', encoding='utf-8') as outf:
            json.dump(out_maps['aaDeletions'], outf)
        CODE

        # gather runtime metrics
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: "3 GB"
        cpu:    4
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x4"
        maxRetries: 2
    }
    output {
        Map[String,String] nextclade_clade   = read_json("NEXTCLADE_CLADE.json")
        Map[String,String] aa_subs_csv       = read_json("NEXTCLADE_AASUBS.json")
        Map[String,String] aa_dels_csv       = read_json("NEXTCLADE_AADELS.json")
        Array[String]      genome_ids        = read_lines("IDLIST")
        String             nextclade_version = read_string("VERSION")
        File               nextalign_msa     = "~{basename}.nextalign.msa.fasta"
        File               nextclade_json    = "~{basename}.nextclade.json"
        File               auspice_json      = "~{basename}.nextclade.auspice.json"
        File               nextclade_tsv     = "~{basename}.nextclade.tsv"
        Int                max_ram_gb        = ceil(read_float("MEM_BYTES")/1000000000)
        Int                runtime_sec       = ceil(read_float("UPTIME_SEC"))
        String             cpu_load          = read_string("CPU_LOAD")
    }
}

task nextmeta_prep {
  input {
    File   gisaid_meta
    File   assembly_meta
    String out_name
    File?  filter_to_ids
  }
  command <<<
    python3 << CODE
    import os.path
    import csv

    # load inputs
    samples = []
    sample_to_gisaid = {}
    sample_to_assembly = {}
    with open('~{gisaid_meta}', 'rt') as inf:
      for row in csv.DictReader(inf):
        s = row['covv_virus_name'][len('hCoV-19/'):]
        samples.append(s)
        sample_to_gisaid[s] = row
    with open('~{assembly_meta}', 'rt') as inf:
      for row in csv.DictReader(inf, delimiter='\t'):
        sample_to_assembly[row['sample']] = row
    filter_to_ids = "~{default='' filter_to_ids}"
    if filter_to_ids:
        with open(filter_to_ids, 'rt') as inf:
            keep_list = set(x.strip() for x in inf)
    else:
        keep_list = None

    # write outputs
    out_headers = ('strain', 'date', 'region', 'country', 'division', 'location', 'length', 'host', 'Nextstrain_clade', 'pango_lineage', 'originating_lab', 'submitting_lab', 'authors', 'purpose_of_sequencing')

    with open('~{out_name}', 'wt') as outf:
      writer = csv.DictWriter(outf, out_headers, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
      writer.writeheader()

      for sample in samples:
        if not filter_to_ids or sample in keep_list:
            geoloc = sample_to_gisaid[sample]['covv_location'].split(' / ')
            writer.writerow({
            'strain': sample,
            'host': 'Human',
            'length': sample_to_assembly[sample]['assembly_length_unambiguous'],
            'Nextstrain_clade': sample_to_assembly[sample]['nextclade_clade'],
            'pango_lineage': sample_to_assembly[sample]['pango_lineage'],
            'region': geoloc[0],
            'country': geoloc[1] if len(geoloc)>1 else '',
            'division': geoloc[2] if len(geoloc)>2 else '',
            'location': geoloc[3] if len(geoloc)>3 else '',
            'date': sample_to_gisaid[sample]['covv_collection_date'],
            'originating_lab': sample_to_gisaid[sample]['covv_orig_lab'],
            'submitting_lab': sample_to_gisaid[sample]['covv_subm_lab'],
            'authors': sample_to_gisaid[sample]['covv_authors'],
            'purpose_of_sequencing': sample_to_gisaid[sample].get('covv_sampling_strategy', sample_to_gisaid[sample].get('covv_add_host_info')),
            })

    CODE
  >>>
  output {
    File nextmeta_tsv = "~{out_name}"
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task derived_cols {
    meta {
        description: "Create derivative columns in nextstrain metadata file (optionally). TSV input and output files may be gzipped."
    }
    input {
        File          metadata_tsv
        String?       lab_highlight_loc
        Array[File]   table_map = []

        String        docker = "quay.io/broadinstitute/viral-core:2.3.1"
        Int           disk_size = 50
    }
    parameter_meta {
        lab_highlight_loc: {
          description: "This option copies the 'originating_lab' and 'submitting_lab' columns to new ones including a prefix, but only if they match certain criteria. The value of this string must be of the form prefix;col_header=value:col_header=value. For example, 'MA;country=USA:division=Massachusetts' will copy the originating_lab and submitting_lab columns to MA_originating_lab and MA_submitting_lab, but only for those rows where country=USA and division=Massachusetts."
        }
        table_map: {
            description: "Mapping tables. Each mapping table is a tsv with a header. The first column is the output column name for this mapping (it will be created or overwritten). The subsequent columns are matching criteria. The value in the first column is written to the output column. The exception is in the case where all match columns are '*' -- in this case, the value in the first column is the column header name to copy over.",
            patterns: ["*.txt", "*.tsv"]
        }
    }
    String basename = basename(basename(metadata_tsv, ".txt"), ".tsv")
    command <<<
        python3<<CODE
        import csv, gzip

        def open_or_gzopen(*args, **kwargs):
            return gzip.open(*args, **kwargs) if args[0].endswith('.gz') else open(*args, **kwargs)

        class Adder_Table_Map:
            def __init__(self, tab_file):
                self._mapper = {}
                self._default_col = None
                with open_or_gzopen(tab_file, 'rt') as inf:
                    reader = csv.DictReader(inf, delimiter='\t')
                    self._col_name = reader.fieldnames[0]
                    self._orig_cols = reader.fieldnames[1:]
                    for row in reader:
                        if all(v=='*' for k,v in row.items() if k in self._orig_cols):
                            self._default_col = row.get(self._col_name)
                        else:
                            k = self._make_key_str(row)
                            v = row.get(self._col_name, '')
                            print("setting {}={} if {}".format(self._col_name, v, k))
                            self._mapper[k] = v
            def _make_key_str(self, row):
                key_str = ':'.join('='.join((k,row.get(k,''))) for k in self._orig_cols)
                return key_str
            def extra_headers(self):
                return (self._col_name,)
            def modify_row(self, row):
                k = self._make_key_str(row)
                v = self._mapper.get(k)
                if v is None and self._default_col:
                   v = row.get(self._default_col, '')
                row[self._col_name] = v
                return row

        class Adder_Source_Lab_Subset:
            def __init__(self, restrict_string):
                self._prefix = restrict_string.split(';')[0]
                self._restrict_map = dict(kv.split('=') for kv in restrict_string.split(';')[1].split(':'))
            def extra_headers(self):
                return (self._prefix + '_originating_lab', self._prefix + '_submitting_lab')
            def modify_row(self, row):
                if all((row.get(k) == v) for k,v in self._restrict_map.items()):
                    row[self._prefix + '_originating_lab'] = row['originating_lab']
                    row[self._prefix + '_submitting_lab']  = row['submitting_lab']
                return row

        def tsv_derived_cols(in_tsv, out_tsv, table_map=None, lab_highlight_loc=None):
            adders = []
            if table_map:
                for t in table_map:
                    adders.append(Adder_Table_Map(t))
            if lab_highlight_loc:
               adders.append(Adder_Source_Lab_Subset(lab_highlight_loc))

            with open_or_gzopen(in_tsv, 'rt') as inf:
                reader = csv.DictReader(inf, delimiter='\t')
                out_headers = reader.fieldnames
                for adder in adders:
                    out_headers.extend(adder.extra_headers())

                with open_or_gzopen(out_tsv, 'wt') as outf:
                    writer = csv.DictWriter(outf, out_headers, delimiter='\t')
                    writer.writeheader()
                    for row in reader:
                        for adder in adders:
                            adder.modify_row(row)
                        writer.writerow(row)

        lab_highlight_loc = "~{default='' lab_highlight_loc}"
        table_map = list(x for x in "~{sep='*' table_map}".split('*') if x)
        tsv_derived_cols(
            "~{metadata_tsv}",
            "~{basename}.derived_cols.txt",
            table_map = table_map,
            lab_highlight_loc = lab_highlight_loc if lab_highlight_loc else None
        )

        CODE
    >>>
    runtime {
        docker: docker
        memory: "1 GB"
        cpu:    1
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File derived_metadata = "~{basename}.derived_cols.txt"
    }
}

task filter_segments {
    input {
        File  all_samples_fasta
        Int   segment = 1
        File? pre_assembled_samples_fasta

        Int   machine_mem_gb = 3
        Int   disk_size = 375
    }
    command <<<
    python3 <<CODE

        segment = "-"+'~{segment}'
        segment_fasta = ""

        with open('~{all_samples_fasta}', 'r') as fasta:
            records=fasta.read().split(">")

            for r in records:
                if len(r.split("\n")) > 1:
                    header = r.split("\n")[0]

                    if segment in header:
                        new_header = header.replace(segment, "")
                        contents = r.replace(header, new_header)
                        segment_fasta += ">"+contents

        if '~{pre_assembled_samples_fasta}':
            with open('~{pre_assembled_samples_fasta}', 'r') as pre_assembled:
                segment_fasta += pre_assembled.read()
        print(segment_fasta)

    CODE
    >>>
    runtime {
        docker: "python:slim"
        memory: machine_mem_gb + " GB"
        cpu:    1
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File assembly_of_segment = stdout()
    }
}

task nextstrain_build_subsample {
    meta {
        description: "Filter and subsample a sequence set using a Nextstrain 'build.yaml' file. See https://docs.nextstrain.org/en/latest/tutorials/SARS-CoV-2/steps/customizing-analysis.html#custom-subsampling-schemes"
    }
    input {
        File   alignment_msa_fasta
        File   sample_metadata_tsv
        String build_name
        File?  builds_yaml
        File?  parameters_yaml
        File?  keep_list
        File?  drop_list

        Int    machine_mem_gb = 50
        String docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        String nextstrain_ncov_repo_commit = "30435fb9ec8de2f045167fb90adfec12f123e80a"
        Int    disk_size = 750
    }
    parameter_meta {
        alignment_msa_fasta: {
          description: "Multiple sequence alignment (aligned fasta)",
          patterns: ["*.fasta"]
        }
        sample_metadata_tsv: {
            description: "Tab-separated metadata file that contain binning variables and values. Must contain all samples: output will be filtered to the IDs present in this file.",
            patterns: ["*.txt", "*.tsv"]
        }
        build_name: {
            description: "Name of the nextstrain 'build' definition to use for subsampling, as defined in the builds_yaml file."
        }
        builds_yaml: {
            description: "YAML-formatted nextstrain 'build' definitions. See https://docs.nextstrain.org/en/latest/tutorials/SARS-CoV-2/steps/customizing-analysis.html#custom-subsampling-schemes for details. If this is not specified, we will use this as a default: https://github.com/nextstrain/ncov/blob/master/my_profiles/example/builds.yaml",
            patterns: ["*.yaml"]
        }
        parameters_yaml: {
            description: "YAML-formatted nextstrain parameter override definitions. If this is not specified, we will use this as a default: https://github.com/nextstrain/ncov/blob/master/defaults/parameters.yaml",
            patterns: ["*.yaml"]
        }
    }
    command <<<
        set -e -o pipefail
        augur version > VERSION

        # pull the ncov repo w/Snakemake rules at pinned version
        wget -q "https://github.com/nextstrain/ncov/archive/~{nextstrain_ncov_repo_commit}.tar.gz"
        mkdir -p ncov
        tar -xf "~{nextstrain_ncov_repo_commit}.tar.gz" -C ncov --strip-components=1
        cd ncov

        # set the config file
        cat > my_profiles/config.yaml <<CONFIG
        configfile:
          - ~{default="defaults/parameters.yaml" parameters_yaml}
          - builds.yaml
        printshellcmds: True
        show-failed-logs: True
        reason: True
        stats: stats.json
        CONFIG

        # point to input data in builds.yaml
        cat > builds.yaml <<INPUTS
        inputs:
            - name: dataset
              sequences: "~{alignment_msa_fasta}"
              metadata: "~{sample_metadata_tsv}"
        INPUTS
        if [ -n "~{default='' builds_yaml}" ]; then
            # user specifies own builds, merge with pointers to data
            cat "~{default='' builds_yaml}" >> builds.yaml
        else
            # user does not specify their own builds.yaml file, so use example file from ncov repo
            # strip away default pointers to example data in example builds
            tail +20 my_profiles/example/builds.yaml >> builds.yaml
        fi

        # hard inclusion list
        # This prepares the "defaults/include.txt" file expected by the subsample Snakemake rule 
        # operating with the default configuration
        KEEP_LIST="~{default='' keep_list}"
        if [ -n "$KEEP_LIST" ]; then
            for i in $(cat defaults/include.txt); do echo $i; done > include-ncov-default-cleannewlines.txt
            cat include-ncov-default-cleannewlines.txt $KEEP_LIST > defaults/include.txt
        fi

        # hard exclusion list
        # This prepares the "defaults/exclude.txt" file expected by the subsample Snakemake rule 
        # operating with the default configuration
        #
        # See:
        #   https://github.com/nextstrain/ncov/blob/7b8e74d80772641fc656310cd2b83d2f11dde76a/workflow/snakemake_rules/main_workflow.smk#L292
        #   https://github.com/nextstrain/ncov/blob/638470f64b980b60e7bb766866b7faa8b7a7c5aa/defaults/parameters.yaml#L57
        DROP_LIST="~{default='' drop_list}"
        if [ -n "$DROP_LIST" ]; then
            for i in $(cat defaults/exclude.txt); do echo $i; done > exclude-ncov-default-cleannewlines.txt
            cat exclude-ncov-default-cleannewlines.txt $DROP_LIST > defaults/exclude.txt
        fi

        # seed input data (skip some upstream steps in the DAG)
        # strip away anything after a space (esp in fasta headers--they break priorities.py)
        mkdir -p results
        cut -f 1 -d ' ' "~{alignment_msa_fasta}" > results/masked_dataset.fasta
        xz -2 results/masked_dataset.fasta

        # execute snakemake on pre-iqtree target
        RAM_MB=$(cat /proc/meminfo | grep MemTotal | perl -lape 's/MemTotal:\s+(\d+)\d\d\d\s+kB/$1/')
        snakemake \
            -j $(nproc) \
            --resources mem_mb=$RAM_MB \
            --profile my_profiles \
            results/"~{build_name}"/"~{build_name}"_subsampled_sequences.fasta.xz

        # decompress xz's for downstream
        xz -d results/"~{build_name}"/"~{build_name}"_subsampled_sequences.fasta.xz
        xz -d results/"~{build_name}"/"~{build_name}"_subsampled_metadata.tsv.xz

        # grab logs and numbers
        set +o pipefail
        grep . logs/subsample_"~{build_name}"_* > ../augur.filter.logs.txt
        grep \> "results/~{build_name}/~{build_name}_subsampled_sequences.fasta" | wc -l | tee ../OUT_COUNT
        for i in results/"~{build_name}"/sample-*.fasta; do
          group=$(basename $i .fasta | cut -f 2- -d -)
          n=$(grep \> $i | wc -l)
          echo -e "$group"'\t'$n
        done > ../counts_by_group

        # gather runtime metrics
        cd ..
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        cpu :   4
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem3_ssd1_v2_x8"
        maxRetries: 2
    }
    output {
        File            subsampled_msa  = "ncov/results/~{build_name}/~{build_name}_subsampled_sequences.fasta"
        File            subsampled_meta = "ncov/results/~{build_name}/~{build_name}_subsampled_metadata.tsv"
        File            subsample_logs  = "augur.filter.logs.txt"
        File            job_stats_json  = "ncov/stats.json"
        Int             sequences_out   = read_int("OUT_COUNT")
        Map[String,Int] counts_by_group = read_map("counts_by_group")
        String          augur_version   = read_string("VERSION")
        Int             max_ram_gb      = ceil(read_float("MEM_BYTES")/1000000000)
        Int             runtime_sec     = ceil(read_float("UPTIME_SEC"))
        String          cpu_load        = read_string("CPU_LOAD")
    }
}

task nextstrain_ncov_defaults {
    input {
        String nextstrain_ncov_repo_commit = "30435fb9ec8de2f045167fb90adfec12f123e80a"
        String docker                      = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int    disk_size = 50
    }
    command {
        set -e
        wget -q "https://github.com/nextstrain/ncov/archive/~{nextstrain_ncov_repo_commit}.tar.gz"
        tar -xf "~{nextstrain_ncov_repo_commit}.tar.gz" --strip-components=1
    }
    runtime {
        docker: docker
        memory: "1 GB"
        cpu:   1
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File clades_tsv      = "defaults/clades.tsv"
        File lat_longs_tsv   = "defaults/lat_longs.tsv"
        File reference_fasta = "defaults/reference_seq.fasta"
        File reference_gb    = "defaults/reference_seq.gb"
        File ids_include     = "defaults/include.txt"
        File ids_exclude     = "defaults/exclude.txt"
        File auspice_config  = "defaults/auspice_config.json"
    }
}

task nextstrain_deduplicate_sequences {
    meta {
        description: "This uses the Nextstrain sanitize_sequences.py script to deduplicate sequences by sequence ID. If the sequences themselves differ, and error is optionally raised. See: https://github.com/nextstrain/ncov/blob/c4747c1f53cd84baaeacdbd044390604d1af2cfc/scripts/sanitize_sequences.py"
    }

    input {
        File sequences_fasta
        Boolean error_on_seq_diff = false

        String nextstrain_ncov_repo_commit = "30435fb9ec8de2f045167fb90adfec12f123e80a"
        String docker                      = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int disk_size = 750
    }

    parameter_meta {
        sequences_fasta: {
          description: "FASTA file with multiple sequences",
          patterns: ["*.fasta","*.fasta.xz","*.fasta.gz"]
        }
        error_on_seq_diff: {
            description: "If error_on_seq_diff=true, an error will be raised if duplicate sequence IDs exist but have different sequences. By default, use the first occurrence of each duplicated sequence."
        }
    }

    String out_basename = basename(basename(basename(basename(sequences_fasta, '.xz'), '.gz'), '.tar'), '.fasta')
    String out_filename = "~{out_basename}_sequences_deduplicated.fasta"
    command {
        set -e
        ncov_path_prefix="/nextstrain/ncov"
        wget -q "https://github.com/nextstrain/ncov/archive/~{nextstrain_ncov_repo_commit}.tar.gz"
        mkdir -p "$ncov_path_prefix"
        tar -xf "~{nextstrain_ncov_repo_commit}.tar.gz" --strip-components=1 -C "$ncov_path_prefix"

        python3 "$ncov_path_prefix/scripts/sanitize_sequences.py" \
        --sequences "~{sequences_fasta}" \
        ${true="--error-on-duplicate-strains" false="" error_on_seq_diff} \
        --output "~{out_filename}"
    }
    runtime {
        docker: docker
        memory: "7 GB"
        cpu:   1
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File sequences_deduplicated_fasta = "~{out_filename}"
    }
}

task nextstrain_ncov_sanitize_gisaid_data {
    meta {
        description: "Sanitize data downloaded from GISAID for use in Nextstrain/augur. See: https://nextstrain.github.io/ncov/data-prep#curate-data-from-the-full-gisaid-database"
    }

    input {
        File sequences_gisaid_fasta
        File metadata_gisaid_tsv

        String? prefix_to_strip

        String nextstrain_ncov_repo_commit = "30435fb9ec8de2f045167fb90adfec12f123e80a"
        String docker                      = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int    disk_size = 750
    }

    parameter_meta {
        sequences_gisaid_fasta: {
          description: "Multiple sequences downloaded from GISAID",
          patterns: ["*.fasta","*.fasta.xz","*.fasta.gz"]
        }
        metadata_gisaid_tsv: {
            description: "Tab-separated metadata file for sequences downloaded from GISAID and passed in via sequences_gisaid_fasta.",
            patterns: ["*.txt", "*.tsv","*.tsv.xz","*.tsv.gz"]
        }
        prefix_to_strip: {
            description: "String prefix to strip from sequence IDs in both the input fasta and metadata files"
        }
    }

    String out_basename = basename(basename(basename(basename(sequences_gisaid_fasta, '.xz'), '.gz'), '.tar'), '.fasta')
    command {
        set -e
        ncov_path_prefix="/nextstrain/ncov"
        wget -q "https://github.com/nextstrain/ncov/archive/~{nextstrain_ncov_repo_commit}.tar.gz"
        mkdir -p "$ncov_path_prefix"
        tar -xf "~{nextstrain_ncov_repo_commit}.tar.gz" --strip-components=1 -C "$ncov_path_prefix"

        python3 "$ncov_path_prefix/scripts/sanitize_sequences.py" \
        --sequences "~{sequences_gisaid_fasta}" \
        ~{"--strip-prefixes=" + prefix_to_strip} \
        --output "~{out_basename}_sequences_sanitized_for_nextstrain.fasta.gz"

        python3 "$ncov_path_prefix/scripts/sanitize_metadata.py" \
        --metadata "~{metadata_gisaid_tsv}" \
        --parse-location-field Location \
        --rename-fields 'Virus name=strain' 'Accession ID=gisaid_epi_isl' 'Collection date=date' 'Clade=GISAID_clade' 'Pango lineage=pango_lineage' 'Host=host' 'Type=virus' 'Patient age=age' \
        ~{"--strip-prefixes=" + prefix_to_strip} \
        --output "~{out_basename}_metadata_sanitized_for_nextstrain.tsv.gz"
    }
    runtime {
        docker: docker
        memory: "7 GB"
        cpu:   1
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File sequences_gisaid_sanitized_fasta = "~{out_basename}_sequences_sanitized_for_nextstrain.fasta.gz"
        File metadata_gisaid_sanitized_tsv    = "~{out_basename}_metadata_sanitized_for_nextstrain.tsv.gz"
    }
}

task filter_subsample_sequences {
    meta {
        description: "Filter and subsample a sequence set. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/filter.html"
    }
    input {
        File           sequences_fasta
        File           sample_metadata_tsv

        Int?           sequences_per_group
        String?        group_by
        File?          include
        File?          exclude

        Boolean        non_nucleotide = true

        Float?         min_date
        Float?         max_date
        Int?           min_length
        File?          priority
        Int?           subsample_seed
        Array[String]? exclude_where
        Array[String]? include_where

        String         docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int            disk_size = 750
    }
    parameter_meta {
        sequences_fasta: {
          description: "Set of sequences (unaligned fasta or aligned fasta -- one sequence per genome) or variants (vcf format) to subsample using augur filter.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
        sample_metadata_tsv: {
          description: "Metadata in tab-separated text format. See https://nextstrain-augur.readthedocs.io/en/stable/faq/metadata.html for details.",
          patterns: ["*.txt", "*.tsv"]
        }
    }
    String out_fname = sub(sub(basename(sequences_fasta), ".vcf", ".filtered.vcf"), ".fasta$", ".filtered.fasta")
    command <<<
        set -e
        augur version > VERSION

        touch wherefile
        VALS="~{write_lines(select_first([exclude_where, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--exclude-where" >> wherefile
            cat $VALS >> wherefile
        fi
        VALS="~{write_lines(select_first([include_where, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--include-where" >> wherefile
            cat $VALS >> wherefile
        fi
        echo --sequences >> wherefile
        echo "~{sequences_fasta}" >> wherefile

        set -o pipefail
        cat wherefile | grep . | tr '\n' '\0' | xargs -0 -t augur filter \
            --metadata "~{sample_metadata_tsv}" \
            ~{"--min-date " + min_date} \
            ~{"--max-date " + max_date} \
            ~{"--min-length " + min_length} \
            ~{true="--non-nucleotide " false=""  non_nucleotide} \
            ~{"--exclude " + exclude} \
            ~{"--include " + include} \
            ~{"--priority " + priority} \
            ~{"--sequences-per-group " + sequences_per_group} \
            ~{"--group-by " + group_by} \
            ~{"--subsample-seed " + subsample_seed} \
            --output "~{out_fname}" | tee STDOUT
        set +o pipefail

        #cat ~{sequences_fasta} | grep \> | wc -l > IN_COUNT
        grep "strains were dropped during filtering" STDOUT | cut -f 1 -d ' ' > DROP_COUNT
        grep "strains passed all filters" STDOUT | cut -f 1 -d ' ' > OUT_COUNT
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: "15 GB"
        cpu :   4
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x4"
        preemptible: 1
        maxRetries: 2
    }
    output {
        File   filtered_fasta    = out_fname
        String augur_version     = read_string("VERSION")
        Int    sequences_dropped = read_int("DROP_COUNT")
        Int    sequences_out     = read_int("OUT_COUNT")
        Int    max_ram_gb        = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec       = ceil(read_float("UPTIME_SEC"))
        String cpu_load          = read_string("CPU_LOAD")
    }
}

task filter_sequences_to_list {
    meta {
        description: "Filter and subsample a sequence set to a specific list of ids in a text file (one id per line)."
    }
    input {
        File         sequences
        Array[File]? keep_list

        String       out_fname = sub(sub(basename(sequences, ".zst"), ".vcf", ".filtered.vcf"), ".fasta$", ".filtered.fasta")
        # Prior docker image: "nextstrain/base:build-20240209T204939Z"
        String       docker = "quay.io/broadinstitute/viral-core:2.3.1"
        Int          disk_size = 750
    }
    parameter_meta {
        sequences: {
          description: "Set of sequences (unaligned fasta or aligned fasta -- one sequence per genome) or variants (vcf format) to subsample using augur filter.",
          patterns: ["*.fasta", "*.fa", "*.fasta.zst", "*.vcf", "*.vcf.gz"]
        }
        keep_list: {
          description: "List of strain ids.",
          patterns: ["*.txt", "*.tsv"]
        }
    }
    command <<<
        set -e
        KEEP_LISTS="~{sep=' ' select_first([keep_list, []])}"

        if [ -n "$KEEP_LISTS" ]; then
            cat $KEEP_LISTS > keep_list.txt
            if [[ "~{sequences}" = *.vcf || "~{sequences}" = *.vcf.gz ]]; then
                # filter vcf file to keep_list.txt
                echo filtering vcf file
                bcftools view --samples-file keep_list.txt "~{sequences}" -Ou | bcftools filter -i "AC>1" -o "~{out_fname}"
                echo "-1" > DROP_COUNT
            else
                # filter fasta file to keep_list.txt
                echo filtering fasta file
    python3 <<CODE
    import Bio.SeqIO
    import util.file
    keep_list = set()
    with open('keep_list.txt', 'rt') as inf:
        keep_list = set(line.strip() for line in inf)
    n_total = 0
    n_kept = 0
    with util.file.open_or_gzopen('~{sequences}', 'rt') as inf:
        with util.file.open_or_gzopen('~{out_fname}', 'wt') as outf:
            for seq in Bio.SeqIO.parse(inf, 'fasta'):
                n_total += 1
                if seq.id in keep_list:
                    n_kept += 1
                    Bio.SeqIO.write(seq, outf, 'fasta')
    n_dropped = n_total-n_kept
    with open('OUT_COUNT', 'wt') as outf:
        outf.write(str(n_kept)+'\n')
    with open('DROP_COUNT', 'wt') as outf:
        outf.write(str(n_dropped)+'\n')
    CODE
            fi

        else
            cp "~{sequences}" "~{out_fname}"
            echo "0" > DROP_COUNT
        fi

        if [[ "~{sequences}" = *.vcf || "~{sequences}" = *.vcf.gz ]]; then
            bcftools query -l "~{out_fname}" > kept_ids.txt
            bcftools view "~{out_fname}" | grep CHROM | awk -F'\t' '{print NF-9}' > OUT_COUNT
        else
            cat "~{out_fname}" | grep \> | cut -c 2- > kept_ids.txt
            cat kept_ids.txt | wc -l > OUT_COUNT
        fi

        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: "7 GB"
        cpu :   2
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x4"
        preemptible: 1
        maxRetries: 2
    }
    output {
        File   filtered_fasta    = out_fname
        Int    sequences_dropped = read_int("DROP_COUNT")
        Int    sequences_out     = read_int("OUT_COUNT")
        File   ids_kept          = "kept_ids.txt"
        Int    max_ram_gb        = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec       = ceil(read_float("UPTIME_SEC"))
        String cpu_load          = read_string("CPU_LOAD")
    }
}

task mafft_one_chr {
    meta {
        description: "Align multiple sequences from FASTA. Only appropriate for closely related (within 99% nucleotide conservation) genomes. See https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html"
    }
    input {
        File     sequences
        File?    ref_fasta
        String   basename
        Boolean  remove_reference = false
        Boolean  keep_length = true
        Boolean  large = false
        Boolean  memsavetree = false

        String   docker = "quay.io/broadinstitute/viral-phylo:2.1.20.2"
        Int      mem_size = 500
        Int      cpus = 64
        Int      disk_size = 750
    }
    command <<<
        set -e

        # decompress sequences if necessary
        GENOMES="~{sequences}"
        if [[ $GENOMES == *.gz ]]; then
            gzip -dc $GENOMES > uncompressed.fasta
            GENOMES="uncompressed.fasta"
        fi

        touch args.txt

        # boolean options
        echo "~{true='--large' false='' large}" >> args.txt
        echo "~{true='--memsavetree' false='' memsavetree}" >> args.txt
        echo "--auto" >> args.txt
        echo "--preservecase" >> args.txt

        # if ref_fasta is specified, use "closely related" mode
        # see https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html
        if [ -f "~{ref_fasta}" ]; then
            echo --addfragments >> args.txt
            echo "$GENOMES" >> args.txt
            echo "~{ref_fasta}" >> args.txt
        else
            echo "$GENOMES" >> args.txt
        fi

        # mafft align to reference in "closely related" mode
        cat args.txt | grep . | xargs -d '\n' mafft --thread -1 \
            ~{true='--keeplength --mapout' false='' keep_length} \
            > msa.fasta

        REMOVE_REF="~{true='--remove-reference' false='' remove_reference}"
        if [ -n "$REMOVE_REF" -a -f "~{ref_fasta}" ]; then
        # remove reference sequence
        python3 <<CODE
        import Bio.SeqIO
        seq_it = Bio.SeqIO.parse('msa.fasta', 'fasta')
        print("dumping " + str(seq_it.__next__().id))
        Bio.SeqIO.write(seq_it, 'msa_drop_one.fasta', 'fasta')
        CODE
            mv msa_drop_one.fasta "~{basename}_aligned.fasta"
        else
            mv msa.fasta "~{basename}_aligned.fasta"
        fi

        # profiling and stats
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: mem_size + " GB"
        cpu :   cpus
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        preemptible: 0
        dx_instance_type: "mem3_ssd1_v2_x36"
        maxRetries: 2
    }
    output {
        File   aligned_sequences = "~{basename}_aligned.fasta"
        Int    max_ram_gb        = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec       = ceil(read_float("UPTIME_SEC"))
        String cpu_load          = read_string("CPU_LOAD")
    }
}

task mafft_one_chr_chunked {
    meta {
        description: "Align multiple sequences from FASTA. Only appropriate for closely related (within 99% nucleotide conservation) genomes. See https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html"
    }
    input {
        File     sequences
        File?    ref_fasta # if reference is not provided, it is assumed to be the first sequence in the "sequences" input file
        String   basename
        Boolean  remove_reference = false

        Int      batch_chunk_size = 2000
        Int      threads_per_job = 2

        String   docker = "quay.io/broadinstitute/viral-phylo:2.1.20.2"
        Int      mem_size = 32
        Int      cpus = 64
        Int      disk_size = 750
    }
    command <<<
        set -e

        # write out ref
        if [ -f "~{ref_fasta}" ]; then
            # if a reference was specified, copy+use it
            cp "~{ref_fasta}" "~{basename}_ref.fasta"
        else
        # otherwise assume the first sequence in the fasta is the ref and export it
        python3 <<CODE
        from Bio import SeqIO
        with open("~{basename}_ref.fasta", "w") as handle:
            record_iter = SeqIO.parse(open("~{sequences}"), "fasta")
            for record in record_iter:
                SeqIO.write(record, handle, "fasta-2line")
                break
        CODE
        fi

        python3 <<CODE
        import os.path
        from Bio import SeqIO
        from itertools import islice, repeat, starmap, takewhile
        from operator import truth
        def batch_iterator(iterable, n):  # https://stackoverflow.com/a/34935239
            return takewhile(truth, map(tuple, starmap(islice, repeat((iter(iterable), n)))))

        record_iter = SeqIO.parse(open("~{sequences}"), "fasta")
        for i, batch in enumerate(batch_iterator(record_iter, ~{batch_chunk_size})):
            chunk_filename = "sequences_mafft_one_chr_chunked_chunk_%i.fasta" % (i + 1)
            # if we're considering the first sequence and the file and the user
            # did not specify a reference sequence
            if i==0 and not os.path.isfile("~{ref_fasta}"):
                # drop first sequence since it's the ref and we've already saved
                # it separately above
                # then write out the rest in this chunk
                with open(chunk_filename, "w") as handle:
                    for record in batch[1:]:
                        SeqIO.write(record, handle, "fasta-2line")
            else:
                with open(chunk_filename, "w") as handle:
                    SeqIO.write(batch, handle, "fasta-2line")
        CODE


        # GNU Parallel refresher:
        # ",," is the replacement string; values after ":::" are substituted where it appears
        parallel --jobs "$(( $((~{cpus}/~{threads_per_job}))<1 ? 1 : $((~{cpus}/~{threads_per_job})) ))" -I ,, \
          "mafft --6merpair --keeplength --preservecase --thread $(((~{threads_per_job}-1)%~{cpus}+1)) --addfragments ,, ~{basename}_ref.fasta > $(basename ,,).msa_chunk.fasta \
          " \
          ::: $(ls -1 sequences_mafft_one_chr_chunked_chunk_*.fasta)

        python3 <<CODE
        import glob, os.path, string
        import Bio.SeqIO
        with open("~{basename}_aligned.fasta", "w") as handle:
            # write a single reference if we are not removing the reference sequence
            if "~{true='remove-ref' false='' remove_reference}" != "remove-ref":
                ref_seq = Bio.SeqIO.parse("~{basename}_ref.fasta", 'fasta')
                Bio.SeqIO.write(ref_seq, handle, 'fasta-2line')
            sorted_alignment_chunks = sorted(glob.glob('*.msa_chunk.fasta'), key=lambda e: int(e.strip(string.ascii_letters+'_.-')))
            for msa_chunk in sorted_alignment_chunks:
                seq_it = Bio.SeqIO.parse(msa_chunk, 'fasta')
                print("dumping " + str(seq_it.__next__().id)) # iterate to remove the reference from each chunk
                Bio.SeqIO.write(seq_it, handle, 'fasta-2line')
        CODE

        # profiling and stats
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: mem_size + " GB"
        cpu :   cpus
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        preemptible: 0
        dx_instance_type: "mem3_ssd1_v2_x36"
        maxRetries: 2
    }
    output {
        File   aligned_sequences = "~{basename}_aligned.fasta"
        Int    max_ram_gb        = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec       = ceil(read_float("UPTIME_SEC"))
        String cpu_load          = read_string("CPU_LOAD")
    }
}

task augur_mafft_align {
    meta {
        description: "Align multiple sequences from FASTA. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/align.html"
    }
    input {
        File    sequences
        File    ref_fasta
        String  basename

        File?   existing_alignment
        Boolean fill_gaps = true
        Boolean remove_reference = true

        String  docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int     disk_size = 750
        Int     mem_size = 180
        Int     cpus = 64
    }
    command <<<
        set -e
        augur version > VERSION
        augur align --sequences "~{sequences}" \
            --reference-sequence "~{ref_fasta}" \
            --output "~{basename}_aligned.fasta" \
            ~{true="--fill-gaps" false="" fill_gaps} \
            ~{"--existing-alignment " + existing_alignment} \
            ~{true="--remove-reference" false="" remove_reference} \
            --debug \
            --nthreads auto
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: mem_size + " GB"
        cpu :   cpus
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        preemptible: 0
        dx_instance_type: "mem3_ssd2_v2_x32"
        maxRetries: 2
    }
    output {
        File   aligned_sequences = "~{basename}_aligned.fasta"
        Int    max_ram_gb        = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec       = ceil(read_float("UPTIME_SEC"))
        String cpu_load          = read_string("CPU_LOAD")
        String augur_version     = read_string("VERSION")
    }
}

task snp_sites {
    input {
        File    msa_fasta
        Boolean allow_wildcard_bases = true
        String  docker = "quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0"
        Int     disk_size = 750
    }
    String out_basename = basename(msa_fasta, ".fasta")
    command {
        snp-sites -V > VERSION
        snp-sites -v ~{true="" false="-c" allow_wildcard_bases} -o "~{out_basename}.vcf" "~{msa_fasta}"
    }
    runtime {
        docker: docker
        memory: "31 GB"
        cpu :   2
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        preemptible: 0
        dx_instance_type: "mem3_ssd1_v2_x4"
        maxRetries: 2
    }
    output {
        File   snps_vcf          = "~{out_basename}.vcf"
        String snp_sites_version = read_string("VERSION")
    }
}

task augur_mask_sites {
    meta {
        description: "Mask unwanted positions from alignment or SNP table. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/mask.html"
    }
    input {
        File   sequences
        File?  mask_bed

        String docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int    disk_size = 750
    }
    parameter_meta {
        sequences: {
          description: "Set of alignments (fasta format) or variants (vcf format) to mask.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_fname = sub(sub(basename(sequences), ".vcf", ".masked.vcf"), ".fasta$", ".masked.fasta")
    command <<<
        set -e
        augur version > VERSION
        BEDFILE=~{select_first([mask_bed, "/dev/null"])}
        if [ -s "$BEDFILE" ]; then
            augur mask --sequences "~{sequences}" \
                --mask "$BEDFILE" \
                --output "~{out_fname}"
        else
            cp "~{sequences}" "~{out_fname}"
        fi
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: "3 GB"
        cpu :   4
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        preemptible: 1
        dx_instance_type: "mem1_ssd1_v2_x4"
        maxRetries: 2
    }
    output {
        File   masked_sequences = out_fname
        Int    max_ram_gb       = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec      = ceil(read_float("UPTIME_SEC"))
        String cpu_load         = read_string("CPU_LOAD")
        String augur_version    = read_string("VERSION")
    }
}

task draft_augur_tree {
    meta {
        description: "Build a tree using iqTree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/tree.html"
    }
    input {
        File    msa_or_vcf

        String  method = "iqtree"
        String  substitution_model = "GTR"
        File?   exclude_sites
        File?   vcf_reference
        String? tree_builder_args

        Int     cpus = 64
        Int     machine_mem_gb = 32
        String  docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int     disk_size = 1250
    }
    parameter_meta {
        msa_or_vcf: {
          description: "Set of alignments (fasta format) or variants (vcf format) to construct a tree from using augur tree (iqTree).",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_basename = basename(basename(basename(msa_or_vcf, '.gz'), '.vcf'), '.fasta')
    command <<<
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=100000 augur tree --alignment "~{msa_or_vcf}" \
            --output "~{out_basename}_~{method}.nwk" \
            --method "~{method}" \
            --substitution-model ~{default="GTR" substitution_model} \
            ~{"--exclude-sites " + exclude_sites} \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--tree-builder-args " + tree_builder_args} \
            --nthreads auto
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        cpu:    cpus
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x36"
        preemptible: 0
        maxRetries: 2
    }
    output {
        File   aligned_tree  = "~{out_basename}_~{method}.nwk"
        Int    max_ram_gb    = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec   = ceil(read_float("UPTIME_SEC"))
        String cpu_load      = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task refine_augur_tree {
    meta {
        description: "Refine an initial tree using sequence metadata and Treetime. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/refine.html"
    }
    input {
        File     raw_tree
        File     msa_or_vcf
        File     metadata

        Int?     gen_per_year
        Float?   clock_rate
        Float?   clock_std_dev
        Boolean  keep_root = true
        String?  root
        Boolean? covariance
        Boolean  keep_polytomies = false
        Int?     precision
        Boolean  date_confidence = true
        String?  date_inference = "marginal"
        String?  branch_length_inference
        String?  coalescent
        Int?     clock_filter_iqd = 4
        String?  divergence_units = "mutations"
        File?    vcf_reference

        String   docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int      disk_size = 750
        Int      machine_mem_gb = 75
    }
    parameter_meta {
        msa_or_vcf: {
          description: "Set of alignments (fasta format) or variants (vcf format) to use to guide Treetime.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_basename = basename(basename(basename(msa_or_vcf, '.gz'), '.vcf'), '.fasta')
    command <<<
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=100000 augur refine \
            --tree "~{raw_tree}" \
            --alignment "~{msa_or_vcf}" \
            --metadata "~{metadata}" \
            --output-tree "~{out_basename}_timetree.nwk" \
            --output-node-data "~{out_basename}_branch_lengths.json" \
            --timetree \
            ~{"--clock-rate " + clock_rate} \
            ~{"--clock-std-dev " + clock_std_dev} \
            ~{"--coalescent " + coalescent} \
            ~{"--clock-filter-iqd " + clock_filter_iqd} \
            ~{"--gen-per-year " + gen_per_year} \
            ~{"--root " + root} \
            ~{"--precision " + precision} \
            ~{"--date-inference " + date_inference} \
            ~{"--branch-length-inference " + branch_length_inference} \
            ~{"--divergence-units " + divergence_units} \
            ~{true="--covariance" false="--no-covariance" covariance} \
            ~{true="--keep-root" false="" keep_root} \
            ~{true="--keep-polytomies" false="" keep_polytomies} \
            ~{true="--date-confidence" false="" date_confidence} \
            ~{"--vcf-reference " + vcf_reference}
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        cpu :   2
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem3_ssd1_v2_x8"
        preemptible: 0
        maxRetries: 2
    }
    output {
        File   tree_refined   = "~{out_basename}_timetree.nwk"
        File   branch_lengths = "~{out_basename}_branch_lengths.json"
        Int    max_ram_gb     = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec    = ceil(read_float("UPTIME_SEC"))
        String cpu_load       = read_string("CPU_LOAD")
        String augur_version  = read_string("VERSION")
    }
}

task ancestral_traits {
    meta {
        description: "Infer ancestral traits based on a tree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/traits.html"
    }
    input {
        File          tree
        File          metadata
        Array[String] columns

        Boolean       confidence = true
        File?         weights
        Float?        sampling_bias_correction

        Int           machine_mem_gb = 32
        String        docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int           disk_size = 750
    }
    String out_basename = basename(tree, '.nwk')
    command <<<
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=100000 augur traits \
            --tree "~{tree}" \
            --metadata "~{metadata}" \
            --columns ~{sep=" " columns} \
            --output-node-data "~{out_basename}_ancestral_traits.json" \
            ~{"--weights " + weights} \
            ~{"--sampling-bias-correction " + sampling_bias_correction} \
            ~{true="--confidence" false="" confidence}
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        cpu :   4
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem3_ssd1_v2_x4"
        preemptible: 1
        maxRetries: 2
    }
    output {
        File   node_data_json = "~{out_basename}_ancestral_traits.json"
        Int    max_ram_gb     = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec    = ceil(read_float("UPTIME_SEC"))
        String cpu_load       = read_string("CPU_LOAD")
        String augur_version  = read_string("VERSION")
    }
}

task ancestral_tree {
    meta {
        description: "Infer ancestral sequences based on a tree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/ancestral.html"
    }
    input {
        File     tree
        File     msa_or_vcf

        String   inference = "joint"
        Boolean  keep_ambiguous = false
        Boolean  infer_ambiguous = false
        Boolean  keep_overhangs = false
        File?    vcf_reference
        File?    output_vcf

        String   docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int      disk_size = 300
    }
    parameter_meta {
        msa_or_vcf: {
          description: "Set of alignments (fasta format) or variants (vcf format) to use to guide Treetime.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_basename = basename(basename(basename(msa_or_vcf, '.gz'), '.vcf'), '.fasta')
    command <<<
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=100000 augur ancestral \
            --tree "~{tree}" \
            --alignment "~{msa_or_vcf}" \
            --output-node-data "~{out_basename}_nt_muts.json" \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--output-vcf " + output_vcf} \
            --output-sequences "~{out_basename}_ancestral_sequences.fasta" \
            ~{true="--keep-overhangs" false="" keep_overhangs} \
            --inference ~{default="joint" inference} \
            ~{true="--keep-ambiguous" false="" keep_ambiguous} \
            ~{true="--infer-ambiguous" false="" infer_ambiguous}
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: "50 GB"
        cpu :   4
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem3_ssd1_v2_x8"
        preemptible: 0
        maxRetries: 2
    }
    output {
        File   nt_muts_json  = "~{out_basename}_nt_muts.json"
        File   sequences     = "~{out_basename}_ancestral_sequences.fasta"
        Int    max_ram_gb    = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec   = ceil(read_float("UPTIME_SEC"))
        String cpu_load      = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task translate_augur_tree {
    meta {
        description: "export augur files to json suitable for auspice visualization. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/translate.html"
    }
    input {
        File   tree
        File   nt_muts
        File   genbank_gb

        File?  genes
        File?  vcf_reference_output
        File?  vcf_reference

        String docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int    disk_size = 300
    }
    String out_basename = basename(tree, '.nwk')
    command <<<
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=500000 augur translate --tree "~{tree}" \
            --ancestral-sequences "~{nt_muts}" \
            --reference-sequence "~{genbank_gb}" \
            ~{"--vcf-reference-output " + vcf_reference_output} \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--genes " + genes} \
            --output-node-data ~{out_basename}_aa_muts.json
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: "2 GB"
        cpu :   1
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
        maxRetries: 2
    }
    output {
        File   aa_muts_json = "~{out_basename}_aa_muts.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        String augur_version = read_string("VERSION")
    }
}

task tip_frequencies {
    meta {
        description: "Estimating frequencies for tips. See https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/frequencies.html"
    }
    input {
        File     tree
        File     metadata

        String   method = "kde"

        Float?   min_date
        Float?   max_date
        Int?     pivot_interval
        String?  pivot_interval_units
        Float?   narrow_bandwidth
        Float?   wide_bandwidth
        Float?   proportion_wide
        Float?   minimal_frequency
        Float?   stiffness
        Float?   inertia
        Boolean  censored = false
        Boolean  include_internal_nodes = false

        Int      machine_mem_gb = 64
        String   docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        String   out_basename = basename(tree, '.nwk')
        Int      disk_size = 200
    }
    command <<<
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=100000 augur frequencies \
            --method "~{method}" \
            --tree "~{tree}" \
            --metadata "~{metadata}" \
            ~{'--min-date ' + min_date} \
            ~{'--max-date ' + max_date} \
            ~{'--pivot-interval ' + pivot_interval} \
            ~{'--pivot-interval-units ' + pivot_interval_units} \
            ~{'--narrow-bandwidth ' + narrow_bandwidth} \
            ~{'--wide-bandwidth ' + wide_bandwidth} \
            ~{'--proportion-wide ' + proportion_wide} \
            ~{'--narrow-bandwidth ' + narrow_bandwidth} \
            ~{'--minimal-frequency ' + minimal_frequency} \
            ~{'--stiffness ' + stiffness} \
            ~{'--inertia ' + inertia} \
            ~{true='--censored' false='' censored} \
            ~{true='--include-internal-nodes' false='' include_internal_nodes} \
            --output "~{out_basename}_tip-frequencies.json"

        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        cpu :   4
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem3_ssd2_x4"
        preemptible: 1
        maxRetries: 2
    }
    output {
        File   node_data_json = "~{out_basename}_tip-frequencies.json"
        Int    max_ram_gb     = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec    = ceil(read_float("UPTIME_SEC"))
        String cpu_load       = read_string("CPU_LOAD")
        String augur_version  = read_string("VERSION")
    }
}

task assign_clades_to_nodes {
    meta {
        description: "Assign taxonomic clades to tree nodes based on mutation information"
    }
    input {
        File tree_nwk
        File nt_muts_json
        File? aa_muts_json
        File ref_fasta
        File clades_tsv

        String docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int    disk_size = 300
    }
    String out_basename = basename(basename(tree_nwk, ".nwk"), "_timetree")
    command <<<
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=100000 augur clades \
        --tree "~{tree_nwk}" \
        --mutations "~{nt_muts_json}" ~{'"' + aa_muts_json + '"'} \
        --reference "~{ref_fasta}" \
        --clades "~{clades_tsv}" \
        --output-node-data "~{out_basename}_clades.json"
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: "2 GB"
        cpu :   1
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
        maxRetries: 2
    }
    output {
        File   node_clade_data_json = "~{out_basename}_clades.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        String augur_version      = read_string("VERSION")
    }
}

task augur_import_beast {
    meta {
        description: "Import BEAST tree into files ready for augur export, including a Newick-formatted tree and node data in json format. See both https://nextstrain-augur.readthedocs.io/en/stable/faq/import-beast.html and https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/import.html for details."
    }
    input {
        File    beast_mcc_tree

        Float?  most_recent_tip_date
        String? tip_date_regex
        String? tip_date_format
        String? tip_date_delimiter

        Int     machine_mem_gb = 3
        String  docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int     disk_size = 150
    }
    String tree_basename = basename(beast_mcc_tree, ".tree")
    command <<<
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=100000 augur import beast \
            --mcc "~{beast_mcc_tree}" \
            --output-tree "~{tree_basename}.nwk" \
            --output-node-data "~{tree_basename}.json" \
            ~{"--most-recent-tip-date " + most_recent_tip_date} \
            ~{"--tip-date-regex " + tip_date_regex} \
            ~{"--tip-date-format " + tip_date_format} \
            ~{"--tip-date-delimeter " + tip_date_delimiter}
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        cpu :   2
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
        maxRetries: 2
    }
    output {
        File   tree_newick    = "~{tree_basename}.nwk"
        File   node_data_json = "~{tree_basename}.json"
        Int    max_ram_gb     = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec    = ceil(read_float("UPTIME_SEC"))
        String cpu_load       = read_string("CPU_LOAD")
        String augur_version  = read_string("VERSION")
    }
}

task export_auspice_json {
    meta {
        description: "export augur files to json suitable for auspice visualization. The metadata tsv input is generally required unless the node_data_jsons comprehensively capture all of it. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/export.html"
    }
    input {
        File           auspice_config
        File?          sample_metadata
        File           tree
        Array[File]    node_data_jsons

        File?          lat_longs_tsv
        File?          colors_tsv
        Array[String]? geo_resolutions
        Array[String]? color_by_metadata
        File?          description_md
        Array[String]? maintainers
        String?        title
        Boolean        include_root_sequence = true

        String out_basename = basename(basename(tree, ".nwk"), "_timetree")

        Int    machine_mem_gb = 64
        String docker = "docker.io/nextstrain/base:build-20240209T204939Z"
        Int    disk_size = 300
    }
    
    command <<<
        set -e -o pipefail
        augur version > VERSION
        touch exportargs

        # --node-data
        if [ -n "~{sep=' ' node_data_jsons}" ]; then
            echo "--node-data" >> exportargs
            cat "~{write_lines(node_data_jsons)}" >> exportargs
        fi

        # --geo-resolutions
        VALS="~{write_lines(select_first([geo_resolutions, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--geo-resolutions" >> exportargs;
        fi
        cat $VALS >> exportargs

        # --color-by-metadata
        VALS="~{write_lines(select_first([color_by_metadata, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--color-by-metadata" >> exportargs;
        fi
        cat $VALS >> exportargs

        # --title
        if [ -n "~{title}" ]; then
            echo "--title" >> exportargs
            echo "~{title}" >> exportargs
        fi

        # --maintainers
        VALS="~{write_lines(select_first([maintainers, []]))}"
        if [ -n "$(cat $VALS)" ]; then
            echo "--maintainers" >> exportargs;
        fi
        cat $VALS >> exportargs

        # some mandatory args to ensure the file is never empty
        echo --tree >> exportargs
        echo "~{tree}" >> exportargs
        echo --auspice-config >> exportargs
        echo "~{auspice_config}" >> exportargs

        (export AUGUR_RECURSION_LIMIT=100000; cat exportargs | grep . | tr '\n' '\0' | xargs -0 -t augur export v2 \
            ~{"--metadata " + sample_metadata} \
            ~{"--lat-longs " + lat_longs_tsv} \
            ~{"--colors " + colors_tsv} \
            ~{"--description " + description_md} \
            ~{true="--include-root-sequence " false=""  include_root_sequence} \
            --output "~{out_basename}_auspice.json")
        touch "~{out_basename}_auspice_root-sequence.json"
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        set +o pipefail
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        cpu :   4
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem3_ssd1_v2_x8"
        preemptible: 0
        maxRetries: 2
    }
    output {
        File   virus_json         = "~{out_basename}_auspice.json"
        File   root_sequence_json = "~{out_basename}_auspice_root-sequence.json"
        Int    max_ram_gb         = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec        = ceil(read_float("UPTIME_SEC"))
        String cpu_load           = read_string("CPU_LOAD")
        String augur_version      = read_string("VERSION")
    }
}
