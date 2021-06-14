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
        cat ~{sep=" " infiles} > "~{output_name}"
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
            for infname in "~{sep='*' infiles}".split('*'):
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

task md5sum {
  input {
    File in_file
  }
  command {
    md5sum ${in_file} | cut -f 1 | tee MD5
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
    String         placeholder_for_missing = ""
  }

  command <<<
    python3<<CODE
    import collections
    import csv
    import gzip

    out_basename = '~{out_basename}'
    join_id = '~{id_col}'
    in_tsvs = '~{sep="*" input_tsvs}'.split('*')
    readers = list(
      csv.DictReader(gzip.open(fn, 'rt') if fn.endswith('.gz') else open(fn, 'rt'), delimiter='\t')
      for fn in in_tsvs)

    # prep the output header
    header = []
    for reader in readers:
        header.extend(reader.fieldnames)
    header = list(collections.OrderedDict(((h,0) for h in header)).keys())
    if not join_id or join_id not in header:
        raise Exception()

    # merge everything in-memory
    out_ids = []
    out_row_by_id = {}
    for reader in readers:
        for row in reader:
            row_id = row[join_id]
            row_out = out_row_by_id.get(row_id, {})
            for h in header:
                # prefer non-empty values from earlier files in in_tsvs, populate from subsequent files only if missing
                if row_out.get(h, '~{placeholder_for_missing}') == '~{placeholder_for_missing}':
                    row_out[h] = row.get(h, '~{placeholder_for_missing}')
            out_row_by_id[row_id] = row_out
            out_ids.append(row_id)
    out_ids = list(collections.OrderedDict(((i,0) for i in out_ids)).keys())

    # write output
    with open(out_basename+'.txt', 'w', newline='') as outf:
        writer = csv.DictWriter(outf, header, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(out_row_by_id[row_id] for row_id in out_ids)
    CODE
  >>>

  output {
    File out_tsv = "${out_basename}.txt"
  }

  runtime {
    memory: "7 GB"
    cpu: 2
    docker: "python:slim"
    disks: "local-disk 100 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2"
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
  }
}

task tsv_stack {
  input {
    Array[File]+ input_tsvs
    String       out_basename
    String       docker = "quay.io/broadinstitute/viral-core:2.1.31"
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
  }
}

task s3_copy {
  input {
    Array[File] infiles
    String      s3_uri_prefix
    File        aws_credentials
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
    cpu: 2
    disks: "local-disk 1000 HDD"
  }
}
