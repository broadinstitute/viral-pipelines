version 1.0

task concatenate {
    meta {
        description: "This is nothing more than unix cat."
    }
    input {
        Array[File] infiles
        String      output_name
        Int         cpus = 4
    }
    command {
        cat ~{sep=" " infiles} > "~{output_name}"
    }
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    cpus
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File combined = "${output_name}"
    }
}

task zcat {
    meta {
        description: "Glue together a bunch of text files that may or may not be compressed (autodetect among gz,xz,bz2,lz4,zst or uncompressed inputs). Optionally compress the output (depending on requested file extension)"
    }
    input {
        Array[File] infiles
        String      output_name
        Int         cpus = 4
    }
    command <<<
        set -e
        python3 <<CODE
        import os.path
        import gzip, lzma, bz2
        import lz4.frame # pypi library: lz4
        import zstandard as zstd # pypi library: zstandard
        import util.file # viral-core

        # magic bytes from here:
        # https://en.wikipedia.org/wiki/List_of_file_signatures
        magic_bytes_to_compressor = {
            b"\x1f\x8b\x08":             gzip.open,      # .gz
            b"\xfd\x37\x7a\x58\x5a\x00": lzma.open,      # .xz
            b"\x42\x5a\x68":             bz2.open,       # .bz2
            b"\x04\x22\x4d\x18":         lz4.frame.open, # .lz4
            b"\x28\xb5\x2f\xfd":         util.file.zstd_open   # .zst (open using function above rather than library function)
        }
        extension_to_compressor = {
            ".gz":   gzip.open,      # .gz
            ".gzip": gzip.open,      # .gz
            ".xz":   lzma.open,      # .xz
            ".bz2":  bz2.open,       # .bz2
            ".lz4":  lz4.frame.open, # .lz4
            ".zst":  util.file.zstd_open,  # .zst (open using function above rather than library function)
            ".zstd": util.file.zstd_open   # .zst (open using function above rather than library function)
        }

        # max number of bytes we need to identify one of the files listed above
        max_len = max(len(x) for x in magic_bytes_to_compressor.keys())

        def open_or_compressed_open(*args, **kwargs):
            input_file = args[0]

            # if the file exists, try to guess the (de) compressor based on "magic numbers"
            # at the very start of the file
            if os.path.isfile(input_file):
                with open(input_file, "rb") as f:
                    file_start = f.read(max_len)
                for magic, compressor_open_fn in magic_bytes_to_compressor.items():
                    if file_start.startswith(magic):
                        print("opening via {}: {}".format(compressor_open_fn.__module__,input_file))
                        return compressor_open_fn(*args, **kwargs)
                # fall back to generic open if compression type could not be determine from magic numbers
                return open(*args, **kwargs)
            else:
                # if this is a new file, try to choose the opener based on file extension
                for ext,compressor_open_fn in extension_to_compressor.items():
                    if str(input_file).lower().endswith(ext):
                        print("opening via {}: {}".format(compressor_open_fn.__module__,input_file))
                        return compressor_open_fn(*args, **kwargs)
                # fall back to generic open if compression type could not be determine from magic numbers
                return open(*args, **kwargs)

        with open_or_compressed_open("~{output_name}", 'wt') as outf:
            for infname in "~{sep='*' infiles}".split('*'):
                with open_or_compressed_open(infname, 'rt') as inf:
                    for line in inf:
                        outf.write(line)
        CODE

        # gather runtime metrics
        cat /proc/uptime | cut -f 1 -d ' ' > UPTIME_SEC
        cat /proc/loadavg > CPU_LOAD
        { cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes || echo 0; } > MEM_BYTES
    >>>
    runtime {
        docker: "quay.io/broadinstitute/viral-core:2.1.33"
        memory: "1 GB"
        cpu:    cpus
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File    combined     = "${output_name}"
        Int     max_ram_gb   = ceil(read_float("MEM_BYTES")/1000000000)
        Int     runtime_sec  = ceil(read_float("UPTIME_SEC"))
        String  cpu_load     = read_string("CPU_LOAD")
    }
}

task sed {
    meta {
        description: "Replace all occurrences of 'search' with 'replace' using sed."
    }
    input {
        File   infile
        String search
        String replace
        String outfilename = "~{basename(infile)}-rename.txt"
    }
    command {
        sed 's/~{search}/~{replace}/g' "~{infile}" > "~{outfilename}"
    }
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 375 LOCAL"
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File outfile = "~{outfilename}"
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
        maxRetries: 2
    }
    output {
        File ids_txt = "~{basename}.txt"
    }
}

task md5sum {
  input {
    File in_file
  }
  command {
    md5sum ~{in_file} | cut -f 1 -d ' ' | tee MD5
  }
  output {
    String md5 = read_string("MD5")
  }
  runtime {
    docker: "ubuntu"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 100 HDD"
    dx_instance_type: "mem1_ssd2_v2_x2"
    maxRetries: 2
  }
}

task fetch_row_from_tsv {
  input {
    File          tsv
    String        idx_col
    String        idx_val
    Array[String] set_default_keys = []
  }
  command <<<
    python3 << CODE
    import csv, gzip, json
    open_or_gzopen = lambda *args, **kwargs: gzip.open(*args, **kwargs) if args[0].endswith('.gz') else open(*args, **kwargs)
    out_dict = {}
    with open_or_gzopen('~{tsv}', 'rt') as inf:
      for row in csv.DictReader(inf, delimiter='\t'):
        if row.get('~{idx_col}') == '~{idx_val}':
          out_dict = row
          break
    for k in '~{sep="*" set_default_keys}'.split('*'):
      if k and k not in out_dict:
        out_dict[k] = ''
    with open('out.json', 'wt') as outf:
      json.dump(out_dict, outf)
    CODE
  >>>
  output {
    Map[String,String] map = read_json('out.json')
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task fetch_col_from_tsv {
  input {
    File          tsv
    String        col
    Boolean       drop_empty = true
    Boolean       drop_header = true
    String        out_name = "~{basename(basename(tsv, '.txt'), '.tsv')}-~{col}.txt"
  }
  command <<<
    python3 << CODE
    import csv, gzip
    col = "~{col}"
    drop_empty = ~{true="True" false="False" drop_empty}
    drop_header = ~{true="True" false="False" drop_header}
    open_or_gzopen = lambda *args, **kwargs: gzip.open(*args, **kwargs) if args[0].endswith('.gz') else open(*args, **kwargs)
    with open_or_gzopen('~{tsv}', 'rt') as inf:
      with open('~{out_name}', 'wt') as outf:
        if not drop_header:
          outf.write(col+'\n')
        for row in csv.DictReader(inf, delimiter='\t'):
          x = row.get(col, '')
          if x or not drop_empty:
            outf.write(x+'\n')
    CODE
  >>>
  output {
    File  out_txt  = "~{out_name}"
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task tsv_join {
  meta {
      description: "Perform a full left outer join on multiple TSV tables. Each input tsv must have a header row, and each must must contain the value of id_col in its header. Inputs may or may not be gzipped. Unix/Mac/Win line endings are tolerated on input, Unix line endings are emitted as output. Unicode text safe."
  }

  input {
    Array[File]+   input_tsvs
    String         id_col
    String         out_basename = "merged"
    String         out_suffix = ".txt"
    Boolean        prefer_first = true
    Int            machine_mem_gb = 7
  }

  command <<<
    python3<<CODE
    import collections
    import csv
    import os.path
    import gzip, lzma, bz2
    import lz4.frame # pypi library: lz4
    import zstandard as zstd # pypi library: zstandard
    import util.file # viral-core

    # magic bytes from here:
    # https://en.wikipedia.org/wiki/List_of_file_signatures
    magic_bytes_to_compressor = {
        b"\x1f\x8b\x08":             gzip.open,      # .gz
        b"\xfd\x37\x7a\x58\x5a\x00": lzma.open,      # .xz
        b"\x42\x5a\x68":             bz2.open,       # .bz2
        b"\x04\x22\x4d\x18":         lz4.frame.open, # .lz4
        b"\x28\xb5\x2f\xfd":         util.file.zstd_open   # .zst (open using function above rather than library function)
    }
    extension_to_compressor = {
        ".gz":   gzip.open,      # .gz
        ".gzip": gzip.open,      # .gz
        ".xz":   lzma.open,      # .xz
        ".bz2":  bz2.open,       # .bz2
        ".lz4":  lz4.frame.open, # .lz4
        ".zst":  util.file.zstd_open,  # .zst (open using function above rather than library function)
        ".zstd": util.file.zstd_open   # .zst (open using function above rather than library function)
    }

    # max number of bytes we need to identify one of the files listed above
    max_len = max(len(x) for x in magic_bytes_to_compressor.keys())

    def open_or_compressed_open(*args, **kwargs):
        input_file = args[0]

        # if the file exists, try to guess the (de) compressor based on "magic numbers"
        # at the very start of the file
        if os.path.isfile(input_file):
            with open(input_file, "rb") as f:
                file_start = f.read(max_len)
            for magic, compressor_open_fn in magic_bytes_to_compressor.items():
                if file_start.startswith(magic):
                    print("opening via {}: {}".format(compressor_open_fn.__module__,input_file))
                    return compressor_open_fn(*args, **kwargs)
            # fall back to generic open if compression type could not be determine from magic numbers
            return open(*args, **kwargs)
        else:
            # if this is a new file, try to choose the opener based on file extension
            for ext,compressor_open_fn in extension_to_compressor.items():
                if str(input_file).lower().endswith(ext):
                    print("opening via {}: {}".format(compressor_open_fn.__module__,input_file))
                    return compressor_open_fn(*args, **kwargs)
            # fall back to generic open if compression type could not be determine from magic numbers
            return open(*args, **kwargs)

    # prep input readers
    out_basename = '~{out_basename}'
    join_id = '~{id_col}'
    in_tsvs = '~{sep="*" input_tsvs}'.split('*')
    readers = list(
      csv.DictReader(open_or_compressed_open(fn, 'rt'), delimiter='\t')
      for fn in in_tsvs)

    # prep the output header
    header = []
    for reader in readers:
        header.extend(reader.fieldnames)
    header = list(collections.OrderedDict(((h,0) for h in header)).keys())
    if not join_id or join_id not in header:
        raise Exception()

    # merge everything in-memory
    prefer_first = ~{true="True" false="False" prefer_first}
    out_ids = []
    out_row_by_id = {}
    for reader in readers:
        for row in reader:
            row_id = row[join_id]
            row_out = out_row_by_id.get(row_id, {})
            for h in header:
                if prefer_first:
                  # prefer non-empty values from earlier files in in_tsvs, populate from subsequent files only if missing
                  if not row_out.get(h):
                      row_out[h] = row.get(h, '')
                else:
                  # prefer non-empty values from later files in in_tsvs
                  if row.get(h):
                      row_out[h] = row[h]
            out_row_by_id[row_id] = row_out
            out_ids.append(row_id)
    out_ids = list(collections.OrderedDict(((i,0) for i in out_ids)).keys())

    # write output
    with open_or_compressed_open(out_basename+'~{out_suffix}', 'w', newline='') as outf:
        writer = csv.DictWriter(outf, header, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(out_row_by_id[row_id] for row_id in out_ids)
    CODE
  >>>

  output {
    File out_tsv = "~{out_basename}~{out_suffix}"
  }

  runtime {
    memory: "~{machine_mem_gb} GB"
    cpu: 4
    docker: "quay.io/broadinstitute/viral-core:2.1.33"
    disks: "local-disk 100 HDD"
    dx_instance_type: "mem1_ssd1_v2_x4"
    maxRetries: 2
  }
}

task tsv_to_csv {
  input {
    File   tsv
    String out_basename = basename(basename(tsv, '.tsv'), '.txt')
  }

  command <<<
    python3<<CODE
    import csv
    out_basename = '~{out_basename}'
    with open('~{tsv}', 'rt') as inf:
      reader = csv.DictReader(inf, delimiter='\t')
      with open(out_basename+'.csv', 'w', newline='') as outf:
          writer = csv.DictWriter(outf, reader.fieldnames, dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
          writer.writeheader()
          writer.writerows(reader)
    CODE
  >>>

  output {
    File csv = "${out_basename}.csv"
  }

  runtime {
    memory: "2 GB"
    cpu: 1
    docker: "python:slim"
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}


task tsv_drop_cols {
    meta {
        description: "Remove any private IDs prior to public release."
    }
    input {
        File          in_tsv
        Array[String] drop_cols
        String        out_filename = basename(in_tsv, '.tsv') + ".drop.tsv"
        String        docker = "quay.io/broadinstitute/py3-bio:0.1.2"
    }
    command <<<
        set -e
        python3<<CODE
        import pandas as pd
        in_tsv = "~{in_tsv}"
        df = pd.read_csv(in_tsv, sep='\t', dtype=str).dropna(how='all')
        drop_cols = list(x for x in '~{sep="*" drop_cols}'.split('*') if x)
        if drop_cols:
            df.drop(columns=drop_cols, inplace=True)
        df.to_csv("~{out_filename}", sep='\t', index=False)
        CODE
    >>>
    runtime {
        docker: docker
        memory: "2 GB"
        cpu:    1
        disks: "local-disk 50 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File out_tsv = "~{out_filename}"
    }
}


task tsv_stack {
  input {
    Array[File]+ input_tsvs
    String       out_basename
    String       docker = "quay.io/broadinstitute/viral-core:2.1.33"
  }

  command {
    csvstack -t --filenames \
      ${sep=' ' input_tsvs} \
      | tr , '\t' \
      > ${out_basename}.txt
  }

  output {
    File out_tsv = "${out_basename}.txt"
  }

  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "${docker}"
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task cat_except_headers {
  input {
    Array[File]+ infiles
    String       out_filename
  }

  command <<<
    awk 'FNR>1 || NR==1' \
      ~{sep=' ' infiles} \
      > ~{out_filename}
  >>>

  output {
    File out_tsv = out_filename
  }

  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu"
    disks: "local-disk 50 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}
task make_empty_file {
  input {
    String out_filename
  }
  command {
    touch "~{out_filename}"
  }
  output {
    File out = "~{out_filename}"
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu"
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task rename_file {
  input {
    File   infile
    String out_filename
  }
  command {
    ln -s "~{infile}" "~{out_filename}"
  }
  output {
    File out = "~{out_filename}"
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu"
    disks: "local-disk 100 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task today {
  input {
    String? timezone
  }
  meta {
    volatile: true
  }
  command {
    ~{default='' 'export TZ=' + timezone}
    date +"%Y-%m-%d" > TODAY
  }
  output {
    String date = read_string("TODAY")
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "quay.io/broadinstitute/viral-baseimage:0.1.20"
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task s3_copy {
  input {
    Array[File] infiles
    String      s3_uri_prefix
    File        aws_credentials
    Int         disk_gb = 1000
    Int         cpus = 2
    String?     nop_block # optional ignored input just to allow blocking
  }
  meta {
    description: "aws s3 cp"
  }
  command {
    set -e
    S3_PREFIX=$(echo "~{s3_uri_prefix}" | sed 's|/*$||')
    mkdir -p ~/.aws
    cp ~{aws_credentials} ~/.aws/credentials
    touch OUT_URIS
    for f in ~{sep=' ' infiles}; do
      aws s3 cp $f $S3_PREFIX/
      echo "$S3_PREFIX/$(basename $f)" >> OUT_URIS
    done
  }
  output {
    Array[String] out_uris = read_lines("OUT_URIS")
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-baseimage:0.1.20"
    memory: "2 GB"
    cpu: cpus
    disks: "local-disk ~{disk_gb} SSD"
    maxRetries: 2
  }
}

task filter_sequences_by_length {
    meta {
        description: "Filter sequences in a fasta file to enforce a minimum count of non-N bases."
    }
    input {
        File   sequences_fasta
        Int    min_non_N = 1

        String docker = "quay.io/broadinstitute/viral-core:2.1.33"
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
        disks:  "local-disk 700 HDD"
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File filtered_fasta    = out_fname
        Int  sequences_in      = read_int("IN_COUNT")
        Int  sequences_dropped = read_int("DROP_COUNT")
        Int  sequences_out     = read_int("OUT_COUNT")
    }
}
