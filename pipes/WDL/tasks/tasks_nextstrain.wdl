version 1.0

task concatenate {
    meta {
        description: "This is nothing more than unix cat."
    }
    input {
        Array[File] infiles
        String      output_name
    }
    command {
        cat ~{sep=" " infiles} > "${output_name}"
    }
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File combined = "${output_name}"
    }
}

task gzcat {
    meta {
        description: "Glue together a bunch of text files that may or may not be compressed (autodetect among gz or uncompressed inputs). Optionally compress the output (depending on requested file extension)"
    }
    input {
        Array[File] infiles
        String      output_name
    }
    command <<<
        python3 <<CODE
        import gzip
        open_or_gzopen = lambda *args, **kwargs: gzip.open(*args, **kwargs) if args[0].endswith('.gz') else open(*args, **kwargs)
        with open_or_gzopen("~{output_name}", 'wt') as outf:
            for infname in "~{sep=' ' infiles}".split(' '):
                with open_or_gzopen(infname, 'rt') as inf:
                    for line in inf:
                        outf.write(line)
        CODE
    >>>
    runtime {
        docker: "python:slim"
        memory: "1 GB"
        cpu:    2
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File combined = "${output_name}"
    }
}

task nextmeta_prep {
  input {
    File   gisaid_meta
    File   assembly_meta
    String out_name
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
      for row in csv.DictReader(inf, delimiter='\t'):
        s = row['covv_virus_name'][len('hCoV-19/'):]
        samples.append(s)
        sample_to_gisaid[s] = row
    with open('~{assembly_meta}', 'rt') as inf:
      for row in csv.DictReader(inf, delimiter='\t'):
        sample_to_assembly[row['sample']] = row

    # write outputs
    out_headers = ('strain', 'date', 'region', 'country', 'division', 'location', 'length', 'host', 'Nextstrain_clade', 'pangolin_lineage', 'originating_lab', 'submitting_lab', 'authors', 'purpose_of_sequencing')

    with open('~{out_name}', 'wt') as outf:
      writer = csv.DictWriter(outf, out_headers, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
      writer.writeheader()

      for sample in samples:
        geoloc = sample_to_gisaid[sample]['covv_location'].split(' / ')
        writer.writerow({
          'strain': sample,
          'host': 'Human',
          'length': sample_to_assembly[sample]['assembly_length_unambiguous'],
          'Nextstrain_clade': sample_to_assembly[sample]['nextclade_clade'],
          'pangolin_lineage': sample_to_assembly[sample]['pango_lineage'],
          'region': geoloc[0],
          'country': geoloc[1] if len(geoloc)>1 else '',
          'division': geoloc[2] if len(geoloc)>2 else '',
          'location': geoloc[3] if len(geoloc)>3 else '',
          'date': sample_to_gisaid[sample]['covv_collection_date'],
          'originating_lab': sample_to_gisaid[sample]['covv_orig_lab'],
          'submitting_lab': sample_to_gisaid[sample]['covv_subm_lab'],
          'authors': sample_to_gisaid[sample]['covv_authors'],
          'purpose_of_sequencing': sample_to_gisaid[sample]['covv_add_host_info'],
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
  }
}

task derived_cols {
    meta {
        description: "Create derivative columns in nextstrain metadata file (optionally). TSV input and output files may be gzipped."
    }
    input {
        File          metadata_tsv
        String?       lab_highlight_loc
        Array[File]   table_map=[]

        String        docker="quay.io/broadinstitute/viral-core:2.1.19"
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
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File derived_metadata = "~{basename}.derived_cols.txt"
    }
}

task fasta_to_ids {
    meta {
        description: "Return the headers only from a fasta file"
    }
    input {
        File sequences_fasta
    }
    String basename = basename(sequences_fasta, ".fasta")
    command {
        cat "~{sequences_fasta}" | grep \> | cut -c 2- > "~{basename}.txt"
    }
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File ids_txt = "~{basename}.txt"
    }
}

task filter_segments {
    input {
        File  all_samples_fasta
        Int?  segment = 1
        File? pre_assembled_samples_fasta

        Int? machine_mem_gb
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
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
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
        File     alignment_msa_fasta
        File     sample_metadata_tsv
        String   build_name
        File?    builds_yaml
        File?    parameters_yaml

        Int?     machine_mem_gb
        String   docker = "nextstrain/base:build-20210127T135203Z"
        String   nextstrain_ncov_repo_commit = "5dbca8a45a64e39057c22163f154db981f7ed5c1"
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
          - ~{default="my_profiles/example/builds.yaml" builds_yaml}
        config:
          - sequences=~{alignment_msa_fasta}
          - metadata=~{sample_metadata_tsv}
        printshellcmds: True
        show-failed-logs: True
        reason: True
        stats: stats.json
        CONFIG

        # seed input data (skip some upstream steps in the DAG)
        # strip away anything after a space (esp in fasta headers--they break priorities.py)
        mkdir -p results
        cut -f 1 -d ' ' "~{alignment_msa_fasta}" > results/filtered.fasta

        # execute snakemake on pre-iqtree target
        RAM_MB=$(free -m | head -2 | tail -1 | awk '{print $2}')
        snakemake \
            -j $(nproc) \
            --resources mem_mb=$RAM_MB \
            --profile my_profiles \
            results/"~{build_name}"/subsampled_alignment.fasta

        # grab logs and numbers
        set +o pipefail
        grep . logs/subsample_"~{build_name}"_* > ../augur.filter.logs.txt
        grep \> "results/~{build_name}/subsampled_alignment.fasta" | wc -l | tee ../OUT_COUNT
        for i in results/"~{build_name}"/sample-*.fasta; do
          group=$(basename $i .fasta | cut -f 2- -d -)
          n=$(grep \> $i | wc -l)
          echo -e "$group"'\t'$n
        done > ../counts_by_group

        # gather runtime metrics
        cd ..
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 200]) + " GB" # priorities.py on 700k genomes
        cpu :   4
        disks:  "local-disk 375 HDD"
        dx_instance_type: "mem3_ssd1_v2_x16"
    }
    output {
        File   subsampled_msa = "ncov/results/~{build_name}/subsampled_alignment.fasta"
        File   subsample_logs = "augur.filter.logs.txt"
        File   job_stats_json = "ncov/stats.json"
        Int    sequences_out  = read_int("OUT_COUNT")
        Map[String,Int] counts_by_group = read_map("counts_by_group")
        String augur_version  = read_string("VERSION")
        Int    max_ram_gb     = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec    = ceil(read_float("UPTIME_SEC"))
        String cpu_load       = read_string("CPU_LOAD")
    }
}

task nextstrain_ncov_defaults {
    input {
        String   nextstrain_ncov_repo_commit = "5dbca8a45a64e39057c22163f154db981f7ed5c1"
        String   docker = "nextstrain/base:build-20210127T135203Z"
    }
    command {
        set -e
        wget -q "https://github.com/nextstrain/ncov/archive/~{nextstrain_ncov_repo_commit}.tar.gz"
        tar -xf "~{nextstrain_ncov_repo_commit}.tar.gz" --strip-components=1
        cat defaults/clades.tsv defaults/subclades.tsv > clades-with-subclades.tsv
    }
    runtime {
        docker: docker
        memory: "1 GB"
        cpu :   1
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File clades_tsv      = "clades-with-subclades.tsv"
        File lat_longs_tsv   = "defaults/lat_longs.tsv"
        File reference_fasta = "defaults/reference_seq.fasta"
        File reference_gb    = "defaults/reference_seq.gb"
        File ids_include     = "defaults/include.txt"
        File ids_exclude     = "defaults/exclude.txt"
        File auspice_config  = "defaults/auspice_config.json"
    }
}

task filter_subsample_sequences {
    meta {
        description: "Filter and subsample a sequence set. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/filter.html"
    }
    input {
        File     sequences_fasta
        File     sample_metadata_tsv

        Int?     sequences_per_group
        String?  group_by
        File?    include
        File?    exclude

        Boolean  non_nucleotide=true

        Float?   min_date
        Float?   max_date
        Int?     min_length
        File?    priority
        Int?     subsample_seed
        Array[String]?  exclude_where
        Array[String]?  include_where

        String   docker = "nextstrain/base:build-20210127T135203Z"
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
    command {
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
        grep "sequences were dropped during filtering" STDOUT | cut -f 1 -d ' ' > DROP_COUNT
        grep "sequences have been written out to" STDOUT | cut -f 1 -d ' ' > OUT_COUNT
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "15 GB"
        cpu :   4
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x4"
        preemptible: 1
    }
    output {
        File   filtered_fasta    = out_fname
        String augur_version     = read_string("VERSION")
        Int    sequences_dropped = read_int("DROP_COUNT")
        Int    sequences_out     = read_int("OUT_COUNT")
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
    }
}

task filter_sequences_by_length {
    meta {
        description: "Filter sequences in a fasta file to enforce a minimum count of non-N bases."
    }
    input {
        File    sequences_fasta
        Int     min_non_N = 1

        String  docker="quay.io/broadinstitute/viral-core:2.1.19"
    }
    parameter_meta {
        sequences_fasta: {
          description: "Set of sequences in fasta format",
          patterns: ["*.fasta", "*.fa"]
        }
        min_non_N: {
          description: "Minimum number of called bases (non-N, non-gap, A, T, C, G, and other non-N ambiguity codes accepted)"
        }
    }
    String out_fname = sub(basename(sequences_fasta), ".fasta", ".filtered.fasta")
    command <<<
    python3 <<CODE
    import Bio.SeqIO
    import gzip
    n_total = 0
    n_kept = 0
    open_or_gzopen = lambda *args, **kwargs: gzip.open(*args, **kwargs) if args[0].endswith('.gz') else open(*args, **kwargs)
    with open_or_gzopen('~{sequences_fasta}', 'rt') as inf:
        with open_or_gzopen('~{out_fname}', 'wt') as outf:
            for seq in Bio.SeqIO.parse(inf, 'fasta'):
                n_total += 1
                ungapseq = seq.seq.ungap().upper()
                if (len(ungapseq) - ungapseq.count('N')) >= ~{min_non_N}:
                    n_kept += 1
                    Bio.SeqIO.write(seq, outf, 'fasta')
    n_dropped = n_total-n_kept
    with open('IN_COUNT', 'wt') as outf:
        outf.write(str(n_total)+'\n')
    with open('OUT_COUNT', 'wt') as outf:
        outf.write(str(n_kept)+'\n')
    with open('DROP_COUNT', 'wt') as outf:
        outf.write(str(n_dropped)+'\n')
    CODE
    >>>
    runtime {
        docker: docker
        memory: "1 GB"
        cpu :   1
        disks:  "local-disk 300 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File   filtered_fasta    = out_fname
        Int    sequences_in      = read_int("IN_COUNT")
        Int    sequences_dropped = read_int("DROP_COUNT")
        Int    sequences_out     = read_int("OUT_COUNT")
    }
}

task filter_sequences_to_list {
    meta {
        description: "Filter and subsample a sequence set to a specific list of ids in a text file (one id per line)."
    }
    input {
        File          sequences
        Array[File]?  keep_list

        String   docker = "nextstrain/base:build-20210127T135203Z"
    }
    parameter_meta {
        sequences: {
          description: "Set of sequences (unaligned fasta or aligned fasta -- one sequence per genome) or variants (vcf format) to subsample using augur filter.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
        keep_list: {
          description: "List of strain ids.",
          patterns: ["*.txt", "*.tsv"]
        }
    }
    String out_fname = sub(sub(basename(sequences), ".vcf", ".filtered.vcf"), ".fasta$", ".filtered.fasta")
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
                bcftools view "~{out_fname}" | grep CHROM | awk -F'\t' '{print NF-9}' > OUT_COUNT
            else
                # filter fasta file to keep_list.txt
                echo filtering fasta file
    python3 <<CODE
    import Bio.SeqIO
    keep_list = set()
    with open('keep_list.txt', 'rt') as inf:
        keep_list = set(line.strip() for line in inf)
    n_total = 0
    n_kept = 0
    with open('~{sequences}', 'rt') as inf:
        with open('~{out_fname}', 'wt') as outf:
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
            echo "-1" > OUT_COUNT
        fi
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    >>>
    runtime {
        docker: docker
        memory: "7 GB"
        cpu :   2
        disks:  "local-disk 200 HDD"
        dx_instance_type: "mem1_ssd1_v2_x4"
        preemptible: 1
    }
    output {
        File   filtered_fasta    = out_fname
        Int    sequences_dropped = read_int("DROP_COUNT")
        Int    sequences_out     = read_int("OUT_COUNT")
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
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

        String   docker = "quay.io/broadinstitute/viral-phylo:2.1.19.1"
        Int      mem_size = 500
        Int      cpus = 64
    }
    command {
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

        # remove reference sequence
        python3 <<CODE
        import Bio.SeqIO
        seq_it = Bio.SeqIO.parse('msa.fasta', 'fasta')
        print("dumping " + str(seq_it.__next__().id))
        Bio.SeqIO.write(seq_it, 'msa_drop_one.fasta', 'fasta')
        CODE
        REMOVE_REF="~{true='--remove-reference' false='' remove_reference}"
        if [ -n "$REMOVE_REF" -a -f "~{ref_fasta}" ]; then
            mv msa_drop_one.fasta "~{basename}_aligned.fasta"
        else
            mv msa.fasta "~{basename}_aligned.fasta"
        fi

        # profiling and stats
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: mem_size + " GB"
        cpu :   cpus
        disks:  "local-disk 375 LOCAL"
        preemptible: 0
        dx_instance_type: "mem3_ssd1_v2_x36"
    }
    output {
        File   aligned_sequences = "~{basename}_aligned.fasta"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
    }
}

task augur_mafft_align {
    meta {
        description: "Align multiple sequences from FASTA. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/align.html"
    }
    input {
        File     sequences
        File     ref_fasta
        String   basename

        File?    existing_alignment
        Boolean  fill_gaps = true
        Boolean  remove_reference = true

        String   docker = "nextstrain/base:build-20210127T135203Z"
    }
    command {
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
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "180 GB"
        cpu :   64
        disks:  "local-disk 750 LOCAL"
        preemptible: 0
        dx_instance_type: "mem3_ssd2_v2_x32"
    }
    output {
        File   aligned_sequences = "~{basename}_aligned.fasta"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task snp_sites {
    input {
        File   msa_fasta
        Boolean allow_wildcard_bases=true
        String docker = "quay.io/biocontainers/snp-sites:2.5.1--hed695b0_0"
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
        disks:  "local-disk 100 HDD"
        preemptible: 0
        dx_instance_type: "mem3_ssd1_v2_x4"
    }
    output {
        File   snps_vcf = "~{out_basename}.vcf"
        String snp_sites_version = read_string("VERSION")
    }
}

task augur_mask_sites {
    meta {
        description: "Mask unwanted positions from alignment or SNP table. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/mask.html"
    }
    input {
        File     sequences
        File?    mask_bed

        String   docker = "nextstrain/base:build-20210127T135203Z"
    }
    parameter_meta {
        sequences: {
          description: "Set of alignments (fasta format) or variants (vcf format) to mask.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_fname = sub(sub(basename(sequences), ".vcf", ".masked.vcf"), ".fasta$", ".masked.fasta")
    command {
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
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu :   4
        disks:  "local-disk 100 HDD"
        preemptible: 1
        dx_instance_type: "mem1_ssd1_v2_x4"
    }
    output {
        File   masked_sequences = out_fname
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version  = read_string("VERSION")
    }
}

task draft_augur_tree {
    meta {
        description: "Build a tree using iqTree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/tree.html"
    }
    input {
        File     msa_or_vcf

        String   method = "iqtree"
        String   substitution_model = "GTR"
        File?    exclude_sites
        File?    vcf_reference
        String?  tree_builder_args

        Int?     cpus
        String   docker = "nextstrain/base:build-20210127T135203Z"
    }
    parameter_meta {
        msa_or_vcf: {
          description: "Set of alignments (fasta format) or variants (vcf format) to construct a tree from using augur tree (iqTree).",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_basename = basename(basename(basename(msa_or_vcf, '.gz'), '.vcf'), '.fasta')
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur tree --alignment "~{msa_or_vcf}" \
            --output "~{out_basename}_~{method}.nwk" \
            --method "~{method}" \
            --substitution-model ~{default="GTR" substitution_model} \
            ~{"--exclude-sites " + exclude_sites} \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--tree-builder-args " + tree_builder_args} \
            --nthreads auto
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "32 GB"
        cpu:    select_first([cpus, 64])
        disks:  "local-disk 750 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x36"
        preemptible: 0
    }
    output {
        File   aligned_tree = "~{out_basename}_~{method}.nwk"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
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

        String   docker = "nextstrain/base:build-20210127T135203Z"
    }
    parameter_meta {
        msa_or_vcf: {
          description: "Set of alignments (fasta format) or variants (vcf format) to use to guide Treetime.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_basename = basename(basename(basename(msa_or_vcf, '.gz'), '.vcf'), '.fasta')
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur refine \
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
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "50 GB"
        cpu :   2
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem3_ssd1_v2_x8"
        preemptible: 0
    }
    output {
        File   tree_refined  = "~{out_basename}_timetree.nwk"
        File   branch_lengths = "~{out_basename}_branch_lengths.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task ancestral_traits {
    meta {
        description: "Infer ancestral traits based on a tree. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/traits.html"
    }
    input {
        File           tree
        File           metadata
        Array[String]  columns

        Boolean        confidence = true
        File?          weights
        Float?         sampling_bias_correction

        String   docker = "nextstrain/base:build-20210127T135203Z"
    }
    String out_basename = basename(tree, '.nwk')
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur traits \
            --tree "~{tree}" \
            --metadata "~{metadata}" \
            --columns ~{sep=" " columns} \
            --output-node-data "~{out_basename}_ancestral_traits.json" \
            ~{"--weights " + weights} \
            ~{"--sampling-bias-correction " + sampling_bias_correction} \
            ~{true="--confidence" false="" confidence}
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu :   2
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
    output {
        File   node_data_json = "~{out_basename}_ancestral_traits.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
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

        String   docker = "nextstrain/base:build-20210127T135203Z"
    }
    parameter_meta {
        msa_or_vcf: {
          description: "Set of alignments (fasta format) or variants (vcf format) to use to guide Treetime.",
          patterns: ["*.fasta", "*.fa", "*.vcf", "*.vcf.gz"]
        }
    }
    String out_basename = basename(basename(basename(msa_or_vcf, '.gz'), '.vcf'), '.fasta')
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur ancestral \
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
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "50 GB"
        cpu :   4
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem3_ssd1_v2_x8"
        preemptible: 0
    }
    output {
        File   nt_muts_json = "~{out_basename}_nt_muts.json"
        File   sequences    = "~{out_basename}_ancestral_sequences.fasta"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
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

        String docker = "nextstrain/base:build-20210127T135203Z"
    }
    String out_basename = basename(tree, '.nwk')
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur translate --tree "~{tree}" \
            --ancestral-sequences "~{nt_muts}" \
            --reference-sequence "~{genbank_gb}" \
            ~{"--vcf-reference-output " + vcf_reference_output} \
            ~{"--vcf-reference " + vcf_reference} \
            ~{"--genes " + genes} \
            --output-node-data ~{out_basename}_aa_muts.json
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "2 GB"
        cpu :   1
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
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

        String   docker = "nextstrain/base:build-20210127T135203Z"
        String   out_basename = basename(tree, '.nwk')
    }
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur frequencies \
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
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu :   2
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
    output {
        File   node_data_json = "~{out_basename}_tip-frequencies.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task assign_clades_to_nodes {
    meta {
        description: "Assign taxonomic clades to tree nodes based on mutation information"
    }
    input {
        File tree_nwk
        File nt_muts_json
        File aa_muts_json
        File ref_fasta
        File clades_tsv

        String docker = "nextstrain/base:build-20210127T135203Z"
    }
    String out_basename = basename(basename(tree_nwk, ".nwk"), "_timetree")
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur clades \
        --tree "~{tree_nwk}" \
        --mutations "~{nt_muts_json}" "~{aa_muts_json}" \
        --reference "~{ref_fasta}" \
        --clades "~{clades_tsv}" \
        --output-node-data "~{out_basename}_clades.json"
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "2 GB"
        cpu :   1
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
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

        Int?    machine_mem_gb
        String  docker = "nextstrain/base:build-20210127T135203Z"
    }
    String tree_basename = basename(beast_mcc_tree, ".tree")
    command {
        set -e
        augur version > VERSION
        AUGUR_RECURSION_LIMIT=10000 augur import beast \
            --mcc "~{beast_mcc_tree}" \
            --output-tree "~{tree_basename}.nwk" \
            --output-node-data "~{tree_basename}.json" \
            ~{"--most-recent-tip-date " + most_recent_tip_date} \
            ~{"--tip-date-regex " + tip_date_regex} \
            ~{"--tip-date-format " + tip_date_format} \
            ~{"--tip-date-delimeter " + tip_date_delimiter}
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: select_first([machine_mem_gb, 3]) + " GB"
        cpu :   2
        disks:  "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 1
    }
    output {
        File   tree_newick    = "~{tree_basename}.nwk"
        File   node_data_json = "~{tree_basename}.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}

task export_auspice_json {
    meta {
        description: "export augur files to json suitable for auspice visualization. The metadata tsv input is generally required unless the node_data_jsons comprehensively capture all of it. See https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/export.html"
    }
    input {
        File        auspice_config
        File?       sample_metadata
        File        tree
        Array[File] node_data_jsons

        File?          lat_longs_tsv
        File?          colors_tsv
        Array[String]? geo_resolutions
        Array[String]? color_by_metadata
        File?          description_md
        Array[String]? maintainers
        String?        title
        Boolean        include_root_sequence = true

        String out_basename = basename(basename(tree, ".nwk"), "_timetree")

        String docker = "nextstrain/base:build-20210127T135203Z"
    }
    
    command {
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

        (export AUGUR_RECURSION_LIMIT=10000; cat exportargs | grep . | tr '\n' '\0' | xargs -0 -t augur export v2 \
            ~{"--metadata " + sample_metadata} \
            ~{"--lat-longs " + lat_longs_tsv} \
            ~{"--colors " + colors_tsv} \
            ~{"--description " + description_md} \
            ~{true="--include-root-sequence " false=""  include_root_sequence} \
            --output "~{out_basename}_auspice.json")
        touch "~{out_basename}_auspice_root-sequence.json"
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes > MEM_BYTES
    }
    runtime {
        docker: docker
        memory: "3 GB"
        cpu :   2
        disks:  "local-disk 100 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        preemptible: 0
    }
    output {
        File   virus_json = "~{out_basename}_auspice.json"
        File   root_sequence_json = "~{out_basename}_auspice_root-sequence.json"
        Int    max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
        Int    runtime_sec = ceil(read_float("UPTIME_SEC"))
        String cpu_load = read_string("CPU_LOAD")
        String augur_version = read_string("VERSION")
    }
}
