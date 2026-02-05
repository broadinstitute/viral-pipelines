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
    Int disk_size = 375
    command <<<
        cat ~{sep=" " infiles} > "~{output_name}"
    >>>
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    cpus
        disks: "local-disk ~{disk_size} HDD"
        disk: "~{disk_size} GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File combined = "~{output_name}"
    }
}

task unpack_archive_to_bucket_path {
    meta {
        description: "Unpack archive(s) to a target location within a Google Storage bucket"
    }
    input {
        # input archive(s)
        Array[File] input_archive_files

        # destination for extracted files
        String  bucket_path_prefix
        String? out_dir_name
        
        # tar options
        Boolean bypass_disk_and_unpack_directly_to_bucket = false
        Int?    archive_wrapper_directories_to_strip
        String  tar_opts = "-v --ignore-zeros --no-ignore-command-error"
        
        # gcloud storage options
        Boolean clobber_existing = false
        String? gcloud_access_token
        String  gcloud_storage_cp_opts = ""

        # execution and resource requirements
        Int    disk_size      = ceil(3.0 * size(input_archive_files[0], "GB")) + 50
        Int    machine_mem_gb = 8
        String docker         = "quay.io/broadinstitute/viral-ngs:3.0.2-core"
    }

    parameter_meta {
        # data I/O inputs
        input_archive_files: {
            description: "List of input archive files to unpack.",
            patterns: ["*.tar", "*.tar.gz", "*.tgz", "*.tar.bz2", "*.tbz2", "*.tar.xz", "*.txz", "*.tar.lz4", "*.tar.zst"]
        }
        bucket_path_prefix: {
            description: "The path prefix to the Google Storage bucket location where the archive contents will be unpacked. This must begin with the bucket name, should start with 'gs://', and can include as many sub-directories as desired. If not provided and the job is running on a GCP-backed Terra instance, the bucket of the associated workspace will be inferred via introspection."
        }
        out_dir_name: {
            description: "Name of the (sub-)directory to unpack the archive contents to within the bucket prefix specified. If not provided, the contents will be unpacked to the bucket prefix."
        }
        
        # tar params
        bypass_disk_and_unpack_directly_to_bucket: {
            description: "(tar) If true, unpack the archive(s) and pipe the contents directly to the gcloud storage upload process, without writing to the disk between extraction and upload. If enabled, minimal disk space will be used beyond storage needed to localize the specified input archive(s), but the task may take significantly longer as each file is uploaded using an independent gcloud storage invocation."
        }
        archive_wrapper_directories_to_strip: {
            description: "(tar) If specified, tar extraction excludes this many top-level directories. (i.e. if all files of a tarball are containined within a top-level subdirectory, and archive_wrapper_directories_to_strip=1, the files files will be extracted without being placed into a corresponding output sub-directory. Equivalent to the parameter '--strip-components' of GNU tar."
        }
        tar_opts: {
            description: "(tar) Options to pass to GNU tar during extraction. By default includes: '-v --ignore-zeros --no-ignore-command-error'"
        }

        # 'gcloud storage cp' params
        clobber_existing: {
            description: "(gcloud storage cp) If true, overwrite files in the target directory of the bucket if they already exist."
        }
        gcloud_access_token: {
            description: "(gcloud storage cp) Access token for the Google Cloud Storage bucket, for account authorized to write to the bucket specified by 'bucket_path_prefix'. If not provided, the gcloud auth configuration of the execution environment will be obtained via 'gcloud auth print-access-token' for the active authenticated user (on Terra, the service worker/'pet' account)."
        }
        gcloud_storage_cp_opts: {
            description: "(gcloud storage cp) Additional options to pass to the 'gcloud storage cp' command at the time of upload."
        }
        

        # execution and resource requirements
        disk_size: {
            description: "Size of the disk to allocate for the task, in GB. Note that if multiple files are provided to 'input_archive_files', and extracted data is written to the disk (bypass_disk_and_unpack_directly_to_bucket=false), the extracted data from one archive will be removed before extracting and uploading data from the next input archive."
        }
        machine_mem_gb: {
            description: "Memory to allocate for the task, in GB."
        }
        docker: {
            description: "Docker image to use for the task. For this task, the image must provide GNU tar and the google-cloud-cli ('gcloud' command)"
        }
    }

    command <<<
        set -e

        # verify gcloud is installed (it should be, if the default docker image is used)
        if ! command -v gcloud &> /dev/null; then
            echo "ERROR: gcloud is not installed; it is required to authenticate to Google Cloud Storage" >&2
            exit 1
        fi

        if ~{if(defined(gcloud_access_token)) then 'true' else 'false'}; then
            # set access token env var expected by gcloud,
            # if provided by the user
            CLOUDSDK_AUTH_ACCESS_TOKEN="~{gcloud_access_token}"
        else
            CLOUDSDK_AUTH_ACCESS_TOKEN="$(gcloud auth print-access-token)"
        fi
        export CLOUDSDK_AUTH_ACCESS_TOKEN

        # check that the gcloud access token is populated
        if [ -z "${CLOUDSDK_AUTH_ACCESS_TOKEN}" ]; then
            echo "ERROR: gcloud access token not found; it must either be provided via the 'gcloud_access_token' input, or made available within the execution environment (via 'gcloud auth print-access-token')" >&2
            exit 1
        fi

        # check whether the bucket path prefix begins with "gs://" and if not, 
        # prepend the 'protocol'; also strip leading or trailing slash if present
        # (for flexibility; this way the user can specify the bucket path prefix with or without the protocol)
        bucket_path_prefix=$(echo "~{bucket_path_prefix}" | sed -e 's|^gs://||' -e 's|/$||' -e 's|^/*||' -e 's|^|gs://|')
        
        # check that, excluding the gs:// 'protocol' prefix, the bucket path prefix is not empty
        if [ -z "${bucket_path_prefix/#gs:\/\//}" ]; then
            echo "ERROR: bucket path prefix is empty" >&2
            exit 1
        fi

        # check whether the user can write to the target bucket
        # by trying  a simple write action, since we cannot rely on
        # the user having the permissions needed to view the IAM policies
        # that determine their (write) access to the bucket 
        if ! echo "write_test" | gcloud storage cp --verbosity error - "${bucket_path_prefix}/.tmp/test-write-access.txt" --quiet; then
            echo "ERROR: user does not have write access to the target bucket: ~{bucket_path_prefix}" >&2
            exit 1
        else
            # clean up the test file if the write test was successful
            gcloud storage rm "${bucket_path_prefix}/.tmp/test-write-access.txt"
        fi

        # for each of the input archives provided, extract the contents to the target bucket
        # either directly via pipe, or from an intermediate location on disk
        for input_archive in ~{sep=' ' input_archive_files}; do
            echo "Processing archive: $(basename "${input_archive}")"

            # if the user has requested to bypass writing to disk between extraction and upload
            if ~{if(bypass_disk_and_unpack_directly_to_bucket) then 'true' else 'false'}; then
                echo "Unpacking archive(s) and piping directly to gcloud storage upload processes (bypassing the disk)..."

                # TODO: parallelize if needed and if the increased memory usage is acceptable
                #       either via GNU parallel ( https://www.gnu.org/software/parallel/parallel_examples.html )
                #       or by simply pushing the tar processes to the background
                
                # pipe each file to a command via stdout, relying GNU tar to pass file information
                # out of band via special environment variables set for each file when using the --to-command
                #
                #   documentation here:
                #     https://www.gnu.org/software/tar/manual/html_section/extract-options.html#Writing-to-an-External-Program
                tar ~{tar_opts} -x \
                    ~{if(defined(archive_wrapper_directories_to_strip)) then "--strip-components=~{archive_wrapper_directories_to_strip}" else ""} \
                    --to-command='gcloud storage cp ~{gcloud_storage_cp_opts} ~{if clobber_existing then "" else "--no-clobber"} --verbosity error - '"${bucket_path_prefix}~{if(defined(out_dir_name)) then '/~{out_dir_name}' else ''}/"'${TAR_REALNAME}' \
                    -f "${input_archive}"
            
            # otherwise extract to disk and then upload to the bucket
            else
                echo 'Extracting archive '"$(basename "${input_archive}")"' to disk before upload...'

                # create a temporary directory to extract the archive contents to
                mkdir -p extracted_tmp

                # extract the archive to the temporary directory
                tar ~{tar_opts} -x \
                --directory "./extracted_tmp" \
                ~{if(defined(archive_wrapper_directories_to_strip)) then "--strip-components=~{archive_wrapper_directories_to_strip}" else ""} \
                -f "${input_archive}"

                pushd extracted_tmp

                echo "Uploading extracted files to the target bucket..."
                
                # gcloud storage rsync the extracted files to the target bucket in the target directory
                gcloud storage rsync \
                    --recursive \
                    ~{if clobber_existing then "" else "--no-clobber"} \
                    --verbosity warning \
                    ~{gcloud_storage_cp_opts} \
                    ./ "${bucket_path_prefix}~{if(defined(out_dir_name)) then '/~{out_dir_name}' else ''}"

                popd 
                rm -r ./extracted_tmp
            fi
        done
    >>>

    runtime {
        docker: docker
        memory: "~{machine_mem_gb} GB"
        cpu:    16
        disks: "local-disk ~{disk_size} LOCAL"
        disk: "~{disk_size} GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x16"
        preemptible: 0
        maxRetries: 1
    }

    output {
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
    Int disk_size = 375
    command <<<
        set -e
        python3 <<CODE
        import os.path
        import gzip, lzma, bz2
        import lz4.frame # pypi library: lz4
        import zstandard # pypi library: zstandard

        # magic bytes from here:
        # https://en.wikipedia.org/wiki/List_of_file_signatures
        magic_bytes_to_compressor = {
            b"\x1f\x8b\x08":             gzip.open,      # .gz
            b"\xfd\x37\x7a\x58\x5a\x00": lzma.open,      # .xz
            b"\x42\x5a\x68":             bz2.open,       # .bz2
            b"\x04\x22\x4d\x18":         lz4.frame.open, # .lz4
            b"\x28\xb5\x2f\xfd":         zstandard.open  # .zst
        }
        extension_to_compressor = {
            ".gz":   gzip.open,      # .gz
            ".gzip": gzip.open,      # .gz
            ".xz":   lzma.open,      # .xz
            ".bz2":  bz2.open,       # .bz2
            ".lz4":  lz4.frame.open, # .lz4
            ".zst":  zstandard.open, # .zst
            ".zstd": zstandard.open  # .zst
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
        { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.peak ]; then cat /sys/fs/cgroup/memory/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
    >>>
    runtime {
        docker: "quay.io/broadinstitute/viral-ngs:3.0.2-core"
        memory: "1 GB"
        cpu:    cpus
        disks: "local-disk ~{disk_size} LOCAL"
        disk: "~{disk_size} GB" # TES
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
    Int disk_size = 375
    command <<<
        sed 's/~{search}/~{replace}/g' "~{infile}" > "~{outfilename}"
    >>>
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk ~{disk_size} LOCAL"
        disk: "~{disk_size} GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
    }
    output {
        File outfile = "~{outfilename}"
    }
}

task tar_extract {
    meta {
        description: "Extract a tar file"
    }
    input {
        File   tar_file
        Int    disk_size = 375
        String tar_opts = "-z"
    }
    command <<<
        mkdir -p unpack
        cd unpack
        tar -xv ~{tar_opts} -f "~{tar_file}"
    >>>
    runtime {
        docker: "quay.io/broadinstitute/viral-ngs:3.0.2-baseimage"
        memory: "2 GB"
        cpu:    1
        disks: "local-disk ~{disk_size} LOCAL"
        disk: "~{disk_size} GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 2
        preemptible: 1
    }
    output {
        Array[File] files = glob("unpack/*")
    }
}

task download_from_url {
    meta {
        description: "Download a file from a URL. This task exists as a workaround until Terra supports this functionality natively (cromwell already does: https://cromwell.readthedocs.io/en/stable/filesystems/HTTP/). http[s] and ftp supported"
        volatile: true
    }
    input {
        String url_to_download

        String? output_filename
        String? additional_wget_opts
        String  request_method="GET"
        Int     request_max_retries=1

        String? md5_hash_expected
        String? md5_hash_expected_file_url
        Boolean save_response_header_to_file = false

        Int     disk_size = 50
    }

    parameter_meta {
      url_to_download: {
        description: "The URL to download; this is passed to wget"
      }
      
      output_filename: {
        description: "The filename to use for the downloaded file. This is optional, though it can be helpful in the event the server does not advise on a filename via the 'Content-Disposition' header."
      }
      additional_wget_opts: {
        description: "Additional options passed to wget as part of the download command."
      }
      request_method: {
        description: "The request method ('GET', 'POST', etc.) passed to wget. Optional (default: 'GET')"
      }
      request_max_retries: {
        description: "The maximum number of (additional) re-tries to attempt in the event of failed download."
      }
      md5_hash_expected: {
        description: "The (binary-mode) md5 hash expected for the file to download. If provided and the value does not match the md5 hash of the downloaded file, the task will fail. mutually exclusive with md5_hash_expected_file_url"
      }
      md5_hash_expected_file_url: {
        description: "The url of a file containing the (binary-mode) md5 hash expected for the file to download. If provided and the value does not match the md5 hash of the downloaded file, the task will fail. mutually exclusive with md5_hash_expected"
      }
      save_response_header_to_file: {
        description: "If save_response_header_to_file=true, http response headers will be saved to an output file. Only applicable for http[s] URLs."
      }
    }

    String download_subdir_local = "downloaded"
    command <<<
        # enforce that only one source of expected md5 hash can be provided
        ~{if defined(md5_hash_expected) && defined(md5_hash_expected_file_url) then 'echo "The inputs \'md5_hash_expected\' and \'md5_hash_expected_file_url\' cannot both be specified; please provide only one."; exit 1;' else ''}

        mkdir -p "~{download_subdir_local}/tmp"
        
        pushd "~{download_subdir_local}"
        
        # ---- download desired file
        pushd "tmp"

        # if a URL-encoded version of the requested download is needed
        #encoded_url=$(python3 -c "import urllib.parse; print urllib.parse.quote('''~{url_to_download}''')")
        
        # get the desired file using wget
        # --content-disposition = use the file name suggested by the server via the Content-Disposition header
        # --trust-server-names = ...and in the event of a redirect, use the value of the final page rather than that of the original url
        # --save-headers = save the headers sent by the HTTP server to the file, preceding the actual contents, with an empty line as the separator.
        wget \
        --read-timeout 3 --waitretry 30 \
        --no-verbose \
        --method ~{request_method} \
        ~{if defined(output_filename) then "--output-document ~{output_filename}" else ""} \
        --tries ~{request_max_retries} \
        --content-disposition --trust-server-names ~{additional_wget_opts} \
        '~{url_to_download}' \
        ~{if save_response_header_to_file then "--save-headers" else ""} || (echo "ERROR: request to ~{request_method} file from URL failed: ~{url_to_download}"; exit 1)

        # ----

        # get the name of the downloaded file
        downloaded_file_name="$(basename "$(ls -1 | head -n1)")"

        if [ ! -f "$downloaded_file_name" ]; then
            echo "Could not locate downloaded file \"$downloaded_file_name\""
            exit 1
        fi
        
        if [ ! -s "$downloaded_file_name" ]; then
            echo "Downloaded file appears empty: \"$downloaded_file_name\""
            exit 1
        fi

        popd # return to downloaded/

        # (only for http(s)) split http response headers from response body
        # since wget stores both in a single file separated by a couple newlines
        if [[ "~{url_to_download}" =~ ^https?:// ]] && ~{if save_response_header_to_file then "true" else "false"}; then
            echo "Saving response headers separately..."
            csplit -f response -s "tmp/${downloaded_file_name}" $'/^\r$/+1' && \
                mv response00 "../${downloaded_file_name}.headers" && \
                mv response01 "${downloaded_file_name}" && \
                rm "tmp/$downloaded_file_name"
        else
            mv "tmp/${downloaded_file_name}" "${downloaded_file_name}"
        fi
        # alternative python implementation to split response headers from body
        #   via https://stackoverflow.com/a/75483099
        #python3 << CODE
        #if ~{if save_response_header_to_file then "True" else "False"}:
        #    with open("tmp/${downloaded_file_name}", "rb") as f_downloaded:
        #        headers, body = f_downloaded.read().split(b"\r\n\r\n", 1)
        #        # write the response header to a file
        #        with open("${downloaded_file_name}.headers", "wb") as f_headers:
        #            f_headers.write(headers)
        #            f_headers.write(b"\r\n")
        #        # save the file body to its final location
        #        with open("${downloaded_file_name}", "wb") as f:
        #            f.write(body)
        #else:
        #    ## if headers are not being saved, move the file to its final destination
        #    import shutil
        #    shutil.move("tmp/${downloaded_file_name}","${downloaded_file_name}")
        #CODE
        
        rm -r "tmp"

        popd # return to job working directory

        check_md5_sum() {
            # $1 =  md5sum expected
            # $2 =  md5sum of downloaded file
            if [[ "$1" != "$2" ]]; then
                echo "ERROR: md5sum of downloaded file ($2) did not match md5sum expected ($1)";
                exit 1
            fi
        }

        md5sum_of_downloaded=$(md5sum --binary "~{download_subdir_local}/${downloaded_file_name}" | cut -f1 -d' ' | tee MD5_SUM_OF_DOWNLOADED_FILE)

        if ~{if defined(md5_hash_expected) then 'true' else 'false'}; then
            md5_hash_expected="~{md5_hash_expected}"
            check_md5_sum "$md5_hash_expected" "$md5sum_of_downloaded"
        fi
        if ~{if defined(md5_hash_expected_file_url) then 'true' else 'false'}; then
            md5_hash_expected="$(curl --silent ~{md5_hash_expected_file_url} | cut -f1 -d' ')"
            check_md5_sum "$md5_hash_expected" "$md5sum_of_downloaded"
        fi

        # report the file size, in bytes
        printf "Downloaded file size (bytes): " && stat --format=%s  "~{download_subdir_local}/${downloaded_file_name}" | tee SIZE_OF_DOWNLOADED_FILE_BYTES
    >>>
    runtime {
        docker: "quay.io/broadinstitute/viral-ngs:3.0.2-baseimage"
        memory: "2 GB"
        cpu:    1
        disks: "local-disk ~{disk_size} LOCAL"
        disk: "~{disk_size} GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
        maxRetries: 0
        preemptible: 1
    }
    output {
        File  downloaded_response_file    = glob("downloaded/*")[0]
        File? downloaded_response_headers = basename(downloaded_response_file) + ".headers"

        Int    file_size_bytes          = read_int("SIZE_OF_DOWNLOADED_FILE_BYTES")
        String md5_sum_of_response_file = read_string("MD5_SUM_OF_DOWNLOADED_FILE")

        File stdout = stdout()
        File stderr = stderr()
    }
}

task sanitize_fasta_headers {
  input {
    File   in_fasta
    String out_filename = "~{basename(in_fasta, '.fasta')}-sanitized.fasta"
  }
  String docker = "quay.io/broadinstitute/py3-bio:0.1.3"
  Int    disk_size = 375
  command <<<
    python3<<CODE
    import re
    import Bio.SeqIO
    with open('~{in_fasta}', 'rt') as inf:
      with open('~{out_filename}', 'wt') as outf:
        for seq in Bio.SeqIO.parse(inf, 'fasta'):
          seq.id = re.sub(r'[^0-9A-Za-z!_-]', '-', seq.id)
          seq.description = seq.id
          seq.name = seq.id
          Bio.SeqIO.write(seq, outf, 'fasta')
    CODE
  >>>
  output {
    File sanitized_fasta = out_filename
  }
  runtime {
    docker: docker
    memory: "2 GB"
    cpu:    2
    disks: "local-disk ~{disk_size} LOCAL"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 1
  }
}

task fasta_to_ids {
    meta {
        description: "Return the headers only from a fasta file"
    }
    input {
        File sequences_fasta
    }
    Int disk_size = 375
    String basename = basename(sequences_fasta, ".fasta")
    command <<<
        cat "~{sequences_fasta}" | grep \> | cut -c 2- > "~{basename}.txt"
    >>>
    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk ~{disk_size} LOCAL"
        disk: "~{disk_size} GB" # TES
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
  Int disk_size = 100
  command <<<
    md5sum ~{in_file} | cut -f 1 -d ' ' | tee MD5
  >>>
  output {
    String md5 = read_string("MD5")
  }
  runtime {
    docker: "ubuntu"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd2_v2_x2"
    maxRetries: 2
  }
}


task json_dict_to_tsv {
  input {
    String      json_data
    String      out_basename = "out"
  }
  File json_file = write_lines([json_data])
  command <<<
    python3 << CODE
    import csv, json
    with open('~{json_file}', 'rt') as inf:
      data = json.load(inf)
      with open('~{out_basename}.tsv', 'wt') as outf:
        writer = csv.DictWriter(outf, fieldnames=data.keys(), delimiter='\t')
        writer.writeheader()
        writer.writerow(data)
    CODE
  >>>
  output {
    File tsv = "~{out_basename}.tsv"
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task fetch_row_from_tsv {
  input {
    File          tsv
    String        idx_col
    String        idx_val
    Array[String] set_default_keys = []
    Array[String] add_header = []
  }
  Int disk_size = 50
  command <<<
    python3 << CODE
    import csv, gzip, json
    open_or_gzopen = lambda *args, **kwargs: gzip.open(*args, **kwargs) if args[0].endswith('.gz') else open(*args, **kwargs)
    out_dict = {}
    fieldnames = "~{sep='*' add_header}".split("*")
    if not fieldnames:
      fieldnames = None
    with open_or_gzopen('~{tsv}', 'rt') as inf:
      for row in csv.DictReader(inf, delimiter='\t', fieldnames=fieldnames):
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
      disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
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
  Int disk_size = 50
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
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
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

  Int disk_size = 50

  command <<<
    python3<<CODE
    import collections
    import csv
    import os.path
    import gzip, lzma, bz2
    import lz4.frame # pypi library: lz4
    import zstandard # pypi library: zstandard

    # magic bytes from here:
    # https://en.wikipedia.org/wiki/List_of_file_signatures
    magic_bytes_to_compressor = {
        b"\x1f\x8b\x08":             gzip.open,      # .gz
        b"\xfd\x37\x7a\x58\x5a\x00": lzma.open,      # .xz
        b"\x42\x5a\x68":             bz2.open,       # .bz2
        b"\x04\x22\x4d\x18":         lz4.frame.open, # .lz4
        b"\x28\xb5\x2f\xfd":         zstandard.open  # .zst
    }
    extension_to_compressor = {
        ".gz":   gzip.open,      # .gz
        ".gzip": gzip.open,      # .gz
        ".xz":   lzma.open,      # .xz
        ".bz2":  bz2.open,       # .bz2
        ".lz4":  lz4.frame.open, # .lz4
        ".zst":  zstandard.open, # .zst
        ".zstd": zstandard.open  # .zst
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
    docker: "quay.io/broadinstitute/viral-ngs:3.0.2-core"
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x4"
    maxRetries: 2
  }
}

task tsv_to_csv {
  input {
    File   tsv
    String out_basename = basename(basename(tsv, '.tsv'), '.txt')
  }

  Int disk_size = 50

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
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
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
        String        docker = "quay.io/broadinstitute/py3-bio:0.1.3"
    }
    Int disk_size = 50
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
        disks: "local-disk ~{disk_size} HDD"
        disk: "~{disk_size} GB" # TES
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
    String       docker = "quay.io/broadinstitute/viral-ngs:3.0.2-core"
  }

  Int disk_size = 50

  command <<<
    csvstack -t --filenames \
      ~{sep=' ' input_tsvs} \
      | tr , '\t' \
      > ~{out_basename}.txt
  >>>

  output {
    File out_tsv = "~{out_basename}.txt"
  }

  runtime {
    memory: "1 GB"
    cpu: 1
    docker: docker
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task cat_except_headers {
  input {
    Array[File]+ infiles
    String       out_filename
  }

  Int disk_size = 50

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
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}
task make_empty_file {
  input {
    String out_filename
  }
  Int disk_size = 10
  command <<<
    touch "~{out_filename}"
  >>>
  output {
    File out = "~{out_filename}"
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu"
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task rename_file {
  input {
    File   infile
    String out_filename
  }
  Int disk_size = 375
  command <<<
    cp "~{infile}" "~{out_filename}"
  >>>
  output {
    File out = "~{out_filename}"
  }
  runtime {
    memory: "2 GB"
    cpu: 2
    docker: "ubuntu"
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task raise {
  input {
    String message = "error!"
  }
  command <<<
    set -e
    echo "~{message}"
    exit 1
  >>>
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu"
    disks:  "local-disk 30 HDD"
    disk: "30 GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task unique_strings {
  input {
    Array[String]  strings
    String         separator=","
  }
  Int disk_size = 50
  command <<<
    cat ~{write_lines(strings)} | sort | uniq > UNIQUE_OUT
    python3<<CODE
    with open('UNIQUE_OUT', 'rt') as inf:
      rows = [line.strip() for line in inf]
    with open('UNIQUE_OUT_JOIN', 'wt') as outf:
      outf.write('~{separator}'.join(rows) + '\n')
    CODE
  >>>
  output {
    Array[String]  sorted_unique = read_lines("UNIQUE_OUT")
    String         sorted_unique_joined = read_string("UNIQUE_OUT_JOIN")
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "python:slim"
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task unique_arrays {
  input {
    Array[Array[String]]  string_arrays
  }
  Int disk_size = 50
  command <<<
    cat ~{write_tsv(string_arrays)} | sort | uniq > UNIQUE_OUT
  >>>
  output {
    Array[Array[String]]  sorted_unique = read_tsv("UNIQUE_OUT")
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "ubuntu"
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}

task today {
  input {
    String? timezone
  }
  Int disk_size = 10
  meta {
    volatile: true
  }
  command <<<
    ~{default='' 'export TZ=' + timezone}
    date +"%Y-%m-%d" > TODAY
  >>>
  output {
    String date = read_string("TODAY")
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "quay.io/broadinstitute/viral-ngs:3.0.2-baseimage"
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB" # TES
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
  command <<<
    set -e
    S3_PREFIX=$(echo "~{s3_uri_prefix}" | sed 's|/*$||')
    mkdir -p ~/.aws
    cp ~{aws_credentials} ~/.aws/credentials
    touch OUT_URIS
    for f in ~{sep=' ' infiles}; do
      aws s3 cp $f $S3_PREFIX/
      echo "$S3_PREFIX/$(basename $f)" >> OUT_URIS
    done
  >>>
  output {
    Array[String] out_uris = read_lines("OUT_URIS")
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs:3.0.2-baseimage"
    memory: "2 GB"
    cpu: cpus
    disks: "local-disk ~{disk_gb} SSD"
    disk: "~{disk_gb} GB" # TES
    maxRetries: 2
  }
}

task string_split {
  meta {
    description: "split a string by a delimiter"
  }
  input {
    String   joined_string
    String   delimiter
  }
  command <<<
    set -e
    python3<<CODE
    with open('TOKENS', 'wt') as outf:
      for token in "~{joined_string}".split("~{delimiter}"):
        outf.write(token + '\n')
    CODE
  >>>
  output {
    Array[String] tokens = read_lines("TOKENS")
  }
  runtime {
    docker: "python:slim"
    memory: "1 GB"
    cpu: 1
    disks: "local-disk 50 SSD"
    disk: "50 GB" # TES
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

        String docker = "quay.io/broadinstitute/viral-ngs:3.0.2-core"
        Int    disk_size = 750
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
                ungapseq = seq.seq.replace("-","").upper()
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
        disks: "local-disk ~{disk_size} HDD"
        disk: "~{disk_size} GB" # TES
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

task pair_files_by_basename {
  input {
    Array[File] files
    String      left_ext
    String      right_ext
  }
  Int disk_gb = 100
  command <<<
    set -e
    cp ~{sep=' ' files} .
  >>>
  output {
    Array[File] left_files  = glob("*.~{left_ext}")
    Array[File] right_files = glob("*.~{right_ext}")
    Array[Pair[File,File]] file_pairs = zip(left_files, right_files)
  }
  runtime {
    memory: "1 GB"
    cpu: 2
    docker: "ubuntu"
    disks:  "local-disk ~{disk_gb} HDD"
    disk: "~{disk_gb} GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 2
  }
}
