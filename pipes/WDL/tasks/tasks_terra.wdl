version 1.1

task gcs_copy {
  input {
    Array[File] infiles
    String      gcs_uri_prefix
    # TO DO: add an input File? for GCP credentials to allow this to work outside of Terra
  }
  meta {
    description: "gcloud storage cp without additional authentication only works on Terra"
  }
  parameter_meta {
    infiles: {
      description: "Input files",
      localization_optional: true,
      stream: true
    }
  }
  command <<<
    set -e
    gcloud storage cp "~{sep='" "' infiles}" ~{gcs_uri_prefix}
  >>>
  output {
    File logs = stdout()
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-ngs:3.0.22-baseimage"
    memory: "1 GB"
    cpu: 1
  }
}

task check_terra_env {
  input {
    String docker = "quay.io/broadinstitute/viral-ngs:3.0.22-baseimage"
  }
  meta {
    description: "task for inspection of backend to determine whether the task is running on Terra and/or GCP"
    volatile: true
  }
  command <<<
    # set -x # echo commands upon execution [commented out to avoid leaking the gcloud auth token]
    set -e # exit on pipe fail

    # create gcloud-related output file
    touch gcloud_config_info.log
    touch google_project_id.txt

    # create Terra-related output files
    touch user_email.txt
    touch workspace_id.txt
    touch workspace_name.txt
    touch workspace_namespace.txt
    touch workspace_bucket_path.txt
    touch input_table_name.txt
    touch input_row_id.txt
    touch method_version.txt
    touch method_source.txt
    touch method_path.txt
    touch top_level_submission_id.txt

    #touch gcp_created_by_attributes.txt
    touch gcp_instance_metadata.json

    # disable the version update alert messages gcloud sometimes emits when executing any command
    gcloud config set component_manager/disable_update_check true

    # write system environment variables to output file
    env | tee -a env_info.log

    echo "false" > RUNNING_ON_GCP_PAPIv2
    echo "false" > RUNNING_ON_GCP_BATCH

    # check if running on GCP
    if curl -s metadata.google.internal -i | grep -E 'Metadata-Flavor:\s+Google'; then 
      echo "Cloud platform appears to be GCP"; 
      echo "true" > RUNNING_ON_GCP

      GCLOUD_OAUTH_BEARER_TOKEN="$(gcloud auth print-access-token)"

      #curl -s -H "Metadata-Flavor: Google" \
      #  "http://metadata.google.internal/computeMetadata/v1/instance/attributes/created-by" | tee gcp_created_by_attributes.txt

      curl -s -H "Metadata-Flavor: Google" \
        "http://metadata.google.internal/computeMetadata/v1/instance/?recursive=true" | tee gcp_instance_metadata.json

      # if BATCH_JOB_UID has a value the job is running on GCP Batch
      # NOTE: PAPIv2 is deprecated and will be removed in the future
      if [[ -n "$BATCH_JOB_UID" ]] || $(jq -rc '.attributes | has("cloudbatch-job-uid")' gcp_instance_metadata.json); then
        echo "Job appears to be running on GCP Batch"
        echo "true"  > RUNNING_ON_GCP_BATCH
      else
        echo "Job appears to be running on GCP via PAPIv2"
        echo "true"  > RUNNING_ON_GCP_PAPIv2
      fi

      # Additional introspection can be performed on GCP by querying the internal metadata server
      #   for details see:
      #     https://cloud.google.com/compute/docs/metadata/predefined-metadata-keys

      # write gcloud env info to output files
      gcloud info | tee -a gcloud_config_info.log
    else 
      echo "NOT running on GCP";
      echo "false" > RUNNING_ON_GCP
    fi

    GOOGLE_PROJECT_ID="$(gcloud config list --format='value(core.project)')"
    echo "$GOOGLE_PROJECT_ID" > google_project_id.txt

    # check whether gcloud project has a "terra-" prefix
    # to determine if running on Terra
    if case ${GOOGLE_PROJECT_ID} in terra-*) ;; *) false;; esac; then
      # (shell-portable regex conditional)
      echo "Job appears to be running on Terra (GCP project ID: ${GOOGLE_PROJECT_ID})"
      echo "true" > RUNNING_ON_TERRA

      # get user e-mail for Terra account via firecloud API
      curl -s -X 'GET' \
        'https://api.firecloud.org/me?userDetailsOnly=true' \
        -H 'accept: application/json' \
        -H "Authorization: Bearer $GCLOUD_OAUTH_BEARER_TOKEN" > user_info.json

        USER_EMAIL="$(jq -cr '.userEmail' user_info.json | tee user_email.txt)"
    else
      echo "NOT running on Terra"
      echo "false" > RUNNING_ON_TERRA
    fi

    if grep --quiet "true" RUNNING_ON_GCP && grep --quiet "true" RUNNING_ON_TERRA; then
      echo "Running on Terra+GCP"

      # === Determine Terra workspace ID and submission ID for the workspace responsible for this job

      # locate the Terra (de)localiztion scripts by running find on one of several known potential locations
      # the location may/does differ when running on GCP via PAPIv2 or via Google batch
      known_possible_terra_script_locations=(
                                              "/cromwell_root"
                                              "/mnt/disks/cromwell_root"
                                            )
      terra_localization_script_dirpath="$(dirname $(realpath $(find "${known_possible_terra_script_locations[@]}" -maxdepth 3 -iname gcs_delocalization.sh -print -quit)))"


      # Scrape various workflow / workspace info from the localization and delocalization scripts.
      #   from: https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/GvsUtils.wdl#L35-L40
      WORKSPACE_ID="$(sed -n -E 's!.*gs://fc-(secure-)?([^\/]+).*!\2!p' ${terra_localization_script_dirpath}/gcs_delocalization.sh | sort -u | tee workspace_id.txt)"
      echo "WORKSPACE_ID:            ${WORKSPACE_ID}"

      # check that workspace ID is a valid UUID
      if ! [[ "$WORKSPACE_ID" =~ ^[0-9a-f]{8}-([0-9a-f]{4}-){3}[0-9a-f]{12}$ ]]; then
        echo "ERROR: WORKSPACE_ID identified by parsing ${terra_localization_script_dirpath}/gcs_delocalization.sh is not a valid UUID"
        exit 1
      fi

      # bucket path prefix
      #BUCKET_PREFIX="$(sed -n -E 's!.*(gs://(fc-(secure-)?[^\/]+)).*!\1!p' /cromwell_root/gcs_delocalization.sh | sort -u | tee bucket_prefix.txt)"
      #echo "BUCKET_PREFIX: ${BUCKET_PREFIX}"

      # top-level submission ID
      TOP_LEVEL_SUBMISSION_ID="$(sed -n -E 's!.*gs://fc-(secure-)?([^\/]+)/submissions/([^\/]+).*!\3!p' ${terra_localization_script_dirpath}/gcs_delocalization.sh | sort -u | tee top_level_submission_id.txt)"
      echo "TOP_LEVEL_SUBMISSION_ID: ${TOP_LEVEL_SUBMISSION_ID}"

      # workflow job ID within submission
      #WORKFLOW_ID="$(sed -n -E 's!.*gs://fc-(secure-)?([^\/]+)/submissions/([^\/]+)/([^\/]+)/([^\/]+).*!\5!p' /cromwell_root/gcs_delocalization.sh | sort -u)"
      
      # other way to obtain Terra project ID, via scraping rather than from gcloud call used above
      #GOOGLE_PROJECT_ID="$(sed -n -E 's!.*(terra-[0-9a-f]+).*# project to use if requester pays$!\1!p' /cromwell_root/gcs_localization.sh | sort -u)"
      # =======================================

      # === request workspace name AND namespace from API, based on bucket path / ID ===
      curl -s -X 'GET' \
        "https://api.firecloud.org/api/workspaces/id/${WORKSPACE_ID}?fields=workspace.name%2Cworkspace.namespace%2Cworkspace.googleProject" \
        -H 'accept: application/json' \
        -H "Authorization: Bearer $GCLOUD_OAUTH_BEARER_TOKEN" > workspace_info.json


      WORKSPACE_NAME="$(jq -cr '.workspace.name | select (.!=null)' workspace_info.json | tee workspace_name.txt)"
      WORKSPACE_NAME_URL_ENCODED="$(jq -rn --arg x "${WORKSPACE_NAME}" '$x|@uri')"
      WORKSPACE_NAMESPACE="$(jq -cr '.workspace.namespace | select (.!=null)' workspace_info.json | tee workspace_namespace.txt)"
      WORKSPACE_BUCKET="$(echo "gs://fc-${WORKSPACE_ID}" | tee workspace_bucket_path.txt)"

      echo "WORKSPACE_NAME:      ${WORKSPACE_NAME}"
      echo "WORKSPACE_NAMESPACE: ${WORKSPACE_NAMESPACE}"
      echo "WORKSPACE_BUCKET:    ${WORKSPACE_BUCKET}"

          # --- less direct way of obtaining workspace info by matching Terra project ID --
          #     preserved here for potential utility in obtaining workspace info for other projects/workspaces
          # get list of workspaces, limiting the output to only the fields we need
          #curl -s -X 'GET' \
          #'https://api.firecloud.org/api/workspaces?fields=workspace.name%2Cworkspace.namespace%2Cworkspace.bucketName%2Cworkspace.googleProject' \
          #-H 'accept: application/json' \
          #-H "Authorization: Bearer $GCLOUD_OAUTH_BEARER_TOKEN" > workspace_list.json

          # extract workspace name
          #WORKSPACE_NAME=$(jq -cr '.[] | select( .workspace.googleProject == "'${GOOGLE_PROJECT_ID}'" ).workspace | .name' workspace_list.json)
          
          # extract workspace namespace
          #WORKSPACE_NAMESPACE=$(jq -cr '.[] | select( .workspace.googleProject == "'${GOOGLE_PROJECT_ID}'" ).workspace | .namespace' workspace_list.json)
          #WORKSPACE_NAME_URL_ENCODED="$(jq -rn --arg x "${WORKSPACE_NAME}" '$x|@uri')"

          # extract workspace bucket
          #WORKSPACE_BUCKET=$(jq -cr '.[] | select( .workspace.googleProject == "'${GOOGLE_PROJECT_ID}'" ).workspace | .bucketName' workspace_list.json)
          # --- end less direct way of obtaining workspace info ---
      # =======================================


      # === obtain info on job submission inputs (table name, row ID) ===
      touch submission_metadata.json
      curl -s 'GET' \
      "https://api.firecloud.org/api/workspaces/${WORKSPACE_NAMESPACE}/${WORKSPACE_NAME_URL_ENCODED}/submissions/${TOP_LEVEL_SUBMISSION_ID}" \
      -H 'accept: application/json' \
      -H "Authorization: Bearer $GCLOUD_OAUTH_BEARER_TOKEN" > submission_metadata.json

      INPUT_TABLE_NAME="$(jq -cr 'if .submissionEntity == null then "" elif (.workflows | length)==1 then .submissionEntity.entityType else [.workflows[].workflowEntity.entityType] | join(",") end' submission_metadata.json  | tee input_table_name.txt)"
      INPUT_ROW_ID="$(jq -cr 'if .submissionEntity == null then "" elif (.workflows | length)==1 then .submissionEntity.entityName else [.workflows[].workflowEntity.entityName] | join(",") end' submission_metadata.json | tee input_row_id.txt)"

      echo "INPUT_TABLE_NAME: $INPUT_TABLE_NAME"
      echo "INPUT_ROW_ID:     $INPUT_ROW_ID"
      # =======================================

      # === obtain info on workflow version (branch/tag) and source (dockstore, etc.) ===
      curl -s 'GET' \
        "https://rawls.dsde-prod.broadinstitute.org/api/workspaces/${WORKSPACE_NAMESPACE}/${WORKSPACE_NAME_URL_ENCODED}/submissions/${TOP_LEVEL_SUBMISSION_ID}/configuration" \
        -H 'accept: application/json' \
        -H "Authorization: Bearer $GCLOUD_OAUTH_BEARER_TOKEN" > workflow_version_info.json

      # .methodConfigVersion corresponds to snapshot of input/output config (or a method version stored in Broad methods repo?)
      #jq -cr .methodConfigVersion workflow_version_info.json
      METHOD_VERSION="$(jq -cr '.methodRepoMethod.methodVersion | select (.!=null)' workflow_version_info.json | tee method_version.txt)"
      METHOD_SOURCE="$(jq -cr '.methodRepoMethod.sourceRepo | select (.!=null)' workflow_version_info.json | tee method_source.txt)"
      METHOD_PATH="$(jq -cr '.methodRepoMethod.methodPath | select (.!=null)' workflow_version_info.json | tee method_path.txt)"

      echo "METHOD_VERSION: $METHOD_VERSION"
      echo "METHOD_SOURCE:  $METHOD_SOURCE"
      echo "METHOD_PATH:    $METHOD_PATH"
      # =======================================
    else 
      echo "Not running on Terra+GCP"
    fi

    # pretty-print environment details to stdout
    # if wraping is desired, add to 'column' command: --output-width 120 --table-wrap 0
    ###### disable stdout messages until fixed
    #echo "=============================================="
    #find . -maxdepth 1 -type f \( -iname 'RUNNING*' -or -iname '*.txt' \) -exec sh -c 'printf "$(basename $1 .txt)\t$(head -n1 $1)\n"' _ {} \; | sort -k1 -d -t $'\t' | column --separator $'\t' --table --table-right 1 --output-separator $'  '
    #echo "=============================================="

    echo -n'' "MEM_BYTES: "; { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } | tee MEM_BYTES
  >>>
  output {
    Boolean is_running_on_terra    = read_boolean("RUNNING_ON_TERRA")

    Boolean is_backed_by_gcp          = read_boolean("RUNNING_ON_GCP")
    Boolean is_running_via_gcp_batch  = read_boolean("RUNNING_ON_GCP_BATCH")
    Boolean is_running_via_gcp_papiv2 = read_boolean("RUNNING_ON_GCP_PAPIv2")

    String google_project_id       = read_string("google_project_id.txt")

    String user_email              = read_string("user_email.txt")

    String workspace_uuid            = read_string("workspace_id.txt")
    String workspace_name          = read_string("workspace_name.txt")
    String workspace_namespace     = read_string("workspace_namespace.txt")
    String workspace_bucket_path   = read_string("workspace_bucket_path.txt")

    String method_version          = read_string("method_version.txt")
    String method_source           = read_string("method_source.txt")
    String method_path             = read_string("method_path.txt")

    #String gcp_created_by_metadata = read_string("gcp_created_by_attributes.txt")
    File   gcp_instance_metadata   = "gcp_instance_metadata.json"

    String input_table_name        = read_string("input_table_name.txt")
    String input_row_id            = read_string("input_row_id.txt")

    String top_level_submission_id = read_string("top_level_submission_id.txt")

    File env_info                  = "env_info.log"
    File gcloud_config_info        = "gcloud_config_info.log"

    Int  max_ram_gb                = ceil(read_float("MEM_BYTES")/1000000000)
  }
  runtime {
    docker: docker
    memory: "1 GB"
    cpu: 1
    maxRetries: 2
  }
}

task upload_reads_assemblies_entities_tsv {
  input {
    String        workspace_name
    String        terra_project
    File          tsv_file
    Array[String] cleaned_reads_unaligned_bams_string
    File          meta_by_filename_json

    String        docker = "schaluvadi/pathogen-genomic-surveillance:api-wdl"
  }
  command <<<
    set -e

    echo ~{sep="," cleaned_reads_unaligned_bams_string} > cleaned_bam_strings.txt

    python3 /projects/cdc-sabeti-covid-19/create_data_tables.py \
        -t "~{tsv_file}" \
        -p "~{terra_project}" \
        -w "~{workspace_name}" \
        -b cleaned_bam_strings.txt \
        -j "~{meta_by_filename_json}" \
        | perl -lape 's/^.*Check your workspace for new (\S+) table.*/$1/' \
        > TABLES_MODIFIED
  >>>
  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    maxRetries: 0
  }
  output {
    Array[String] tables = read_lines('TABLES_MODIFIED')
  }
}

task upload_entities_tsv {
  input {
    String        workspace_name
    String        terra_project
    File          tsv_file

    String        docker = "quay.io/broadinstitute/viral-ngs:3.0.22-baseimage"
  }
  meta {
    volatile: true
  }
  command <<<
    set -e
    python3<<CODE
    import sys
    from firecloud import api as fapi
    response = fapi.upload_entities_tsv(
      '~{terra_project}', '~{workspace_name}', '~{tsv_file}', model="flexible")
    if response.status_code != 200:
        print('ERROR UPLOADING: See full error message:')
        print(response.text)
        sys.exit(1)
    else:
        print("Upload complete. Check your workspace for new table!")
    CODE
  >>>
  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    maxRetries: 0
  }
  output {
    Array[String] stdout = read_lines(stdout())
  }
}

task download_entities_tsv {
  input {
    String  terra_project
    String  workspace_name
    String  table_name
    String  outname = "~{terra_project}-~{workspace_name}-~{table_name}.tsv"
    String? nop_input_string # this does absolutely nothing, except that it allows an optional mechanism for you to block execution of this step upon the completion of another task in your workflow

    String  docker = "quay.io/broadinstitute/viral-ngs:3.0.22-baseimage"
  }

  meta {
    volatile: true
  }

  command <<<
    python3<<CODE
    import csv
    import json
    import collections

    from firecloud import api as fapi

    workspace_project = '~{terra_project}'
    workspace_name = '~{workspace_name}'
    table_name = '~{table_name}'
    out_fname = '~{outname}'
    nop_string = '~{default="" nop_input_string}'

    # load terra table and convert to list of dicts
    # I've found that fapi.get_entities_tsv produces malformed outputs if funky chars are in any of the cells of the table
    table = json.loads(fapi.get_entities(workspace_project, workspace_name, table_name).text)
    headers = collections.OrderedDict()
    rows = []
    headers[table_name + "_id"] = 0
    for row in table:
        outrow = row['attributes']
        for x in outrow.keys():
            headers[x] = 0
            if type(outrow[x]) == dict and set(outrow[x].keys()) == set(('itemsType', 'items')):
                outrow[x] = outrow[x]['items']
        outrow[table_name + "_id"] = row['name']
        rows.append(outrow)

    # dump to tsv
    with open(out_fname, 'w', newline='') as outf:
      writer = csv.DictWriter(outf, headers.keys(), delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
      writer.writeheader()
      writer.writerows(rows)
    CODE
  >>>
  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    maxRetries: 2
  }
  output {
    File tsv_file = '~{outname}'
  }
}

task create_or_update_sample_tables {
  input {
    String flowcell_run_id

    String workspace_namespace
    String workspace_name

    Array[String]  raw_reads_unaligned_bams
    Array[String]  cleaned_reads_unaligned_bams

    File           meta_by_filename_json
    File?          read_counts_raw_json
    File?          read_counts_cleaned_json

    String  sample_table_name  = "sample"
    String  library_table_name = "library"

    String  docker = "quay.io/broadinstitute/viral-ngs:3.0.22-core"
  }

  meta {
    volatile: true
  }

  command <<<
    set -e
    python3<<CODE
    flowcell_data_id  = '~{flowcell_run_id}'
    workspace_project = '~{workspace_namespace}'
    workspace_name    = '~{workspace_name}'

    # import required packages
    import sys
    import collections
    import json
    import csv
    import re
    import pandas as pd
    import numpy as np
    from firecloud import api as fapi

    # sanitize table names to conform to Terra naming requirements (alphanumeric, underscores, dashes only)
    sample_table_name = re.sub(r'[^a-zA-Z0-9_-]', '_', '~{sample_table_name}')
    library_table_name = re.sub(r'[^a-zA-Z0-9_-]', '_', '~{library_table_name}')
    lib_col_name = f"entity:{library_table_name}_id"

    print(workspace_project + "\n" + workspace_name)

    # process read counts if available
    read_counts_raw = {}
    read_counts_cleaned = {}
    if '~{default="" read_counts_raw_json}':
        with open('~{default="" read_counts_raw_json}','rt') as inf:
            read_counts_raw = json.load(inf)
    if '~{default="" read_counts_cleaned_json}':
        with open('~{default="" read_counts_cleaned_json}','rt') as inf:
            read_counts_cleaned = json.load(inf)

    # create tsv to populate library table with raw_bam and cleaned_bam columns
    raw_bams_list               = '~{sep="*" raw_reads_unaligned_bams}'.split('*')
    raw_library_id_list         = [bam.split("/")[-1].replace(".bam", "") for bam in raw_bams_list]
    df_library_table_raw_bams   = pd.DataFrame({lib_col_name : raw_library_id_list, "raw_bam" : raw_bams_list})

    cleaned_bams_list           = '~{sep="*" cleaned_reads_unaligned_bams}'.split('*')
    cleaned_library_id_list     = [bam.split("/")[-1].replace(".bam", "").replace(".cleaned", "") for bam in cleaned_bams_list]
    df_library_table_clean_bams = pd.DataFrame({lib_col_name : cleaned_library_id_list, "cleaned_bam" : cleaned_bams_list})
    cleaned_bam_names           = set(df_library_table_clean_bams[lib_col_name])

    df_library_bams = pd.merge(df_library_table_raw_bams, df_library_table_clean_bams, on=lib_col_name, how="outer")
    library_bams_tsv = flowcell_data_id + "-all_bams.tsv"
    df_library_bams.to_csv(library_bams_tsv, sep="\t", index=False)
    library_bam_names = set(df_library_bams[lib_col_name])
    print("libraries in bams: {}".format(len(library_bam_names)))

    # load library metadata from demux json / samplesheet
    with open('~{meta_by_filename_json}',"r") as meta_fp:
        library_meta_dict = json.load(meta_fp)

    # create tsv to populate library table with metadata from demux json / samplesheet
    # to do: maybe just merge this into df_library_bams instead and make a single tsv output
    library_meta_fname = "library_metadata.tsv"
    with open(library_meta_fname, 'w', newline='') as outf:
      copy_cols = ["sample_original", "spike_in", "control", "batch_lib", "library", "lane", "library_id_per_sample", "library_strategy", "library_source", "library_selection", "design_description"]
      out_header = [lib_col_name, 'flowcell', 'read_count_raw', 'read_count_cleaned'] + copy_cols
      print(f"library_metadata.tsv output header: {out_header}")
      writer = csv.DictWriter(outf, out_header, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
      writer.writeheader()

      out_rows = []
      for library in library_meta_dict.values():
        if library['run'] in library_bam_names:
          out_row = {col: library.get(col, '') for col in copy_cols}
          out_row[lib_col_name] = library['run']
          out_row['flowcell'] = flowcell_data_id
          out_row['read_count_raw'] = read_counts_raw.get(library['run'], '')
          out_row['read_count_cleaned'] = read_counts_cleaned.get(library['run'], '')
          out_rows.append(out_row)
      writer.writerows(out_rows)

    # grab the meta_by_filename values to create new sample->library mappings
    # restrict to libraries/samples that we actually have cleaned bam files for
    sample_to_libraries = {}
    libraries_in_bams = set()
    for library_id, data in library_meta_dict.items():
        sample_id = data['sample']
        sample_to_libraries.setdefault(sample_id, [])
        if library_id in cleaned_bam_names:
            sample_to_libraries[sample_id].append(library_id)
            libraries_in_bams.add(library_id)
        else:
            print (f"missing {library_id} from bam list")
    print("json describes {} libraries from {} unique samples".format(len(library_meta_dict), len(sample_to_libraries)))
    print("json describes {} libraries we have bam files for and {} libraries we will ignore".format(len(libraries_in_bams), len(library_meta_dict) - len(libraries_in_bams)))

    # API call to get existing sample->library mappings  <-- THIS IS THE VOLATILE PART
    # get_entities -> python list of dicts
    def get_entities_to_table(project, workspace, table_name):
        table = json.loads(fapi.get_entities(project, workspace, table_name).text)
        headers = collections.OrderedDict()
        rows = []
        headers[table_name + "_id"] = 0
        for row in table:
            outrow = row['attributes']
            for x in outrow.keys():
                headers[x] = 0
                if type(outrow[x]) == dict and set(outrow[x].keys()) == set(('itemsType', 'items')):
                    outrow[x] = outrow[x]['items']
            outrow[table_name + "_id"] = row['name']
            rows.append(outrow)
        return (headers, rows)
    header, rows = get_entities_to_table(workspace_project, workspace_name, sample_table_name)
    df_sample = pd.DataFrame.from_records(rows, columns=header, index=sample_table_name + "_id")
    print(df_sample.index)

    # create tsv to populate sample table with new sample->library mappings
    def test_non_empty_value(value):
      # this function exists because pandas / numpy arrays don't behave like python lists with regards to coercion to truth values
      # Check for numpy NaN (which coerces to python True!)
      if isinstance(value, float) and pd.isna(value):
        return False
      # Check for numpy arrays (which refuse to coerce to logical values and throw ValueError instead!)
      if isinstance(value, (np.ndarray, pd.Series, pd.DataFrame)):
        return value.size > 0 and value.any()
      # Default to normal python behavior
      return bool(value)

    sample_fname = 'sample_membership.tsv'
    with open(sample_fname, 'wt') as outf:
        outf.write(f'entity:{sample_table_name}_id\tlibraries\n')
        merged_sample_ids = set()
        for sample_id, libraries in sample_to_libraries.items():
            if sample_id in df_sample.index and "libraries" in df_sample.columns and test_non_empty_value(df_sample.libraries[sample_id]):
                # merge in new sample->library mappings with any pre-existing sample->library mappings
                already_associated_libraries = [entity["entityName"] for entity in df_sample.libraries[sample_id] if entity.get("entityName")]
                libraries = list(set(libraries + already_associated_libraries))
                print (f"\tsample {sample_id} pre-exists in Terra table, merging old members {already_associated_libraries} with new members {libraries}")
                merged_sample_ids.add(sample_id)

            outf.write(f'{sample_id}\t{json.dumps([{"entityType":library_table_name,"entityName":library_name} for library_name in libraries])}\n')
    print(f"wrote {len(sample_to_libraries)} samples to {sample_fname} where {len(merged_sample_ids)} samples were already in the Terra table")

    # write everything to the Terra table! -- TO DO: move this to separate task
    for fname in (library_bams_tsv, library_meta_fname, sample_fname):
        response = fapi.upload_entities_tsv(workspace_project, workspace_name, fname, model="flexible")
        if response.status_code != 200:
            print(f'ERROR UPLOADING {fname}: See full error message:')
            print(response.text)
            sys.exit(1)
        else:
            print("Upload complete. Check your workspace for new table!")

    CODE
    { if [ -f /sys/fs/cgroup/memory.peak ]; then cat /sys/fs/cgroup/memory.peak; elif [ -f /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes; else echo "0"; fi } > MEM_BYTES
  >>>
  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    maxRetries: 2
  }
  output {
    File library_metadata_tsv = "library_metadata.tsv"
    File sample_membership_tsv = "sample_membership.tsv"
    File library_bams_tsv = "~{flowcell_run_id}-all_bams.tsv"
    File stdout_log = stdout()
    File stderr_log = stderr()
    Int  max_ram_gb = ceil(read_float("MEM_BYTES")/1000000000)
  }
}

task find_illumina_files_in_directory {
  input {
    String  illumina_dir
    String? fastq_dir
    Int?    lane
    Boolean include_undetermined = false
    String  docker = "quay.io/broadinstitute/viral-ngs:3.0.22-baseimage"
  }
  parameter_meta {
    illumina_dir: {
      description: "GCS bucket path to Illumina run directory (e.g., gs://bucket/path/to/run)",
      category: "required"
    }
    fastq_dir: {
      description: "Override path to fastq files (defaults to {illumina_dir}/fastq)",
      category: "advanced"
    }
    lane: {
      description: "If specified, filter outputs to only include FASTQs from this sequencing lane number. Fails if the requested lane does not exist on this flowcell. FASTQs whose names carry no lane token (DRAGEN no-lane-splitting output) represent all lanes merged, and are always included.",
      category: "advanced"
    }
    include_undetermined: {
      description: "If true, include Undetermined_* FASTQs (unassigned reads) in the outputs. Defaults to false.",
      category: "advanced"
    }
    runinfo_xml: {
      description: "GCS path to RunInfo.xml file from the Illumina run directory",
      category: "output"
    }
    fastqs: {
      description: "All FASTQ files found, after applying the lane and Undetermined filters. FASTQs whose names do not parse are still included here, but are omitted from raw_reads_fastq_pairs.",
      category: "output"
    }
    raw_reads_fastq_pairs: {
      description: "FASTQ files grouped by sample, as pairs (PE) or singles (SE). Array of arrays where each inner array contains 1 (SE) or 2 (PE) file paths.",
      category: "output"
    }
  }
  Int disk_size = 20
  command <<<
    set -e -o pipefail

    # Strip trailing slashes from illumina_dir to avoid double-slash issues
    ILLUMINA_DIR="~{illumina_dir}"
    ILLUMINA_DIR="$(echo "$ILLUMINA_DIR" | sed 's:/*$::')"

    # Find RunInfo.xml - check base level first, then search recursively.
    # NB: buffer the listing to a file instead of piping it into `head`. Under
    # `set -o pipefail`, `head` closing the pipe early makes gcloud take SIGPIPE and
    # the pipeline reports failure even though the listing succeeded. That false
    # negative gets likelier on the recursive search, which can return many lines.
    echo "Searching for RunInfo.xml at: $ILLUMINA_DIR/RunInfo.xml" >&2
    if gcloud storage ls "$ILLUMINA_DIR/RunInfo.xml" > runinfo_ls.txt 2>gcloud_error.txt && [ -s runinfo_ls.txt ]; then
      echo "Found RunInfo.xml at base level" >&2
    else
      echo "Not found at base level, searching recursively..." >&2
      cat gcloud_error.txt >&2
      if gcloud storage ls "$ILLUMINA_DIR/**/RunInfo.xml" > runinfo_ls.txt 2>gcloud_error.txt && [ -s runinfo_ls.txt ]; then
        echo "Found RunInfo.xml via recursive search" >&2
      else
        echo "ERROR: RunInfo.xml not found in $ILLUMINA_DIR" >&2
        echo "gcloud error output:" >&2
        cat gcloud_error.txt >&2
        exit 1
      fi
    fi

    head -1 runinfo_ls.txt > runinfo_path.txt
    RUNINFO_PATH=$(cat runinfo_path.txt)
    echo "Found RunInfo.xml at: $RUNINFO_PATH"

    # Read LaneCount from RunInfo.xml. This is the authoritative source for how many
    # lanes the flowcell actually has, and lets us reject an out-of-range lane request
    # up front rather than letting it surface later as a confusing "0 fastqs found".
    gcloud storage cat "$RUNINFO_PATH" > runinfo.xml 2>/dev/null || touch runinfo.xml
    LANE_COUNT=""
    if [ -s runinfo.xml ]; then
      LANE_COUNT="$(awk 'match($0, /LaneCount="[0-9]+"/) { s = substr($0, RSTART, RLENGTH); gsub(/[^0-9]/, "", s); print s; exit }' runinfo.xml)"
    fi
    echo "$LANE_COUNT" > lane_count.txt
    echo "RunInfo.xml declares LaneCount=${LANE_COUNT:-unknown}" >&2

    # Locate the fastqs. An explicit fastq_dir always wins and is never second-guessed.
    # Otherwise try the conventional $ILLUMINA_DIR/fastq, which is where flat delivery
    # layouts put them, and fall back to a recursive search -- an on-instrument DRAGEN
    # demux buries them several levels down, e.g. Analysis/1/Data/fastq/.
    touch all_fastqs.txt
    if [ -n "~{fastq_dir}" ]; then
      FASTQ_DIR="~{fastq_dir}"
      FASTQ_DIR="$(echo "$FASTQ_DIR" | sed 's:/*$::')"
      echo "Using supplied fastq_dir: $FASTQ_DIR" >&2
      gcloud storage ls "$FASTQ_DIR/*.fastq.gz" 2>/dev/null > all_fastqs.txt || touch all_fastqs.txt
    else
      FASTQ_DIR="$ILLUMINA_DIR/fastq"
      echo "Looking for fastqs at: $FASTQ_DIR" >&2
      gcloud storage ls "$FASTQ_DIR/*.fastq.gz" 2>/dev/null > all_fastqs.txt || touch all_fastqs.txt

      if [ ! -s all_fastqs.txt ]; then
        echo "None found there; searching recursively under $ILLUMINA_DIR ..." >&2
        gcloud storage ls "$ILLUMINA_DIR/**/*.fastq.gz" 2>/dev/null > recursive_fastqs.txt || touch recursive_fastqs.txt

        # A re-demux leaves several Analysis/N/ trees behind. Settle on exactly one
        # directory rather than blending them, preferring the highest-numbered;
        # sort -V orders naturally, so Analysis/10 sorts above Analysis/2.
        sed 's:/[^/]*$::' recursive_fastqs.txt | sort -u -V > fastq_dirs_all.txt

        # First drop any candidate nested inside another candidate. A directory that
        # contains another one is the real delivery; its subdirectory is incidental,
        # and would otherwise win the sort below purely for being deeper.
        awk 'NR==FNR { d[NR] = $0; n = NR; next }
             { keep = 1
               for (i = 1; i <= n; i++)
                 if (d[i] != $0 && index($0, d[i] "/") == 1) { keep = 0; break }
               if (keep) print }' fastq_dirs_all.txt fastq_dirs_all.txt > fastq_dirs.txt
        NUM_DIRS="$(wc -l < fastq_dirs.txt)"

        if [ "$NUM_DIRS" -eq 0 ]; then
          FASTQ_DIR="$ILLUMINA_DIR (searched recursively)"
        else
          if [ "$NUM_DIRS" -gt 1 ]; then
            echo "WARNING: fastqs found under $NUM_DIRS directories:" >&2
            sed 's/^/  /' fastq_dirs.txt >&2
            echo "WARNING: using the highest-numbered; set fastq_dir to choose another" >&2
          fi
          FASTQ_DIR="$(tail -1 fastq_dirs.txt)"
          echo "Selected fastq directory: $FASTQ_DIR" >&2
          # keep only files directly in FASTQ_DIR, not in directories nested below it
          awk -v d="$FASTQ_DIR" 'index($0, d "/") == 1 && index(substr($0, length(d) + 2), "/") == 0' \
            recursive_fastqs.txt > all_fastqs.txt
        fi
      fi
    fi

    echo "$FASTQ_DIR" > fastq_dir_used.txt
    if [ ! -s all_fastqs.txt ]; then
      echo "WARNING: No fastq.gz files found in $FASTQ_DIR" >&2
    fi

    # Parse fastq filenames and group into pairs/singles
    python3 << 'CODE'
    import re
    import json
    import sys

    # Lane filter (None if not specified)
    lane_filter = ~{if defined(lane) then lane else "None"}
    include_undetermined = ~{if include_undetermined then "True" else "False"}

    def read_first_line(path, default=''):
        try:
            with open(path, 'rt') as f:
                return f.readline().strip()
        except FileNotFoundError:
            return default

    fastq_dir_searched = read_first_line('fastq_dir_used.txt', '(unknown)')

    lane_count_raw = read_first_line('lane_count.txt')
    lane_count = int(lane_count_raw) if lane_count_raw.isdigit() else None

    # Validate the *requested* lane against the run's own LaneCount before filtering
    # anything. A lane that does not exist on this flowcell is an invalid request and
    # should say so, rather than silently reducing to an empty result set.
    if lane_filter is not None:
        if lane_count is None:
            print("WARNING: could not read LaneCount from RunInfo.xml; skipping lane range check")
        elif lane_filter < 1 or lane_filter > lane_count:
            print(f"ERROR: lane {lane_filter} was requested, but this run has {lane_count} lane(s) "
                  f"per LaneCount in RunInfo.xml. Valid lanes are 1..{lane_count}.")
            sys.exit(1)

    # Read all fastq paths
    with open('all_fastqs.txt', 'rt') as f:
        fastqs = [line.strip() for line in f if line.strip()]

    # Two patterns, tried in order. A single regex with an optional (?:_L(\d+))? group
    # will not do: the greedy (.+) swallows "_L001" and skips the optional group, so
    # genuinely lane-split runs would silently stop reporting their lane.
    #   lane-split (BCL Convert default): Sample1_S1_L001_R1_001.fastq.gz
    #   no-lane-splitting (DRAGEN):       Sample1_S1_R1_001.fastq.gz
    pattern_lane   = re.compile(r'^(.+)_L(\d+)_R([12])_(\d+)\.fastq\.gz$')
    pattern_nolane = re.compile(r'^(.+)_R([12])_(\d+)\.fastq\.gz$')

    groups = {}
    unparsed = []
    skipped_undetermined = []
    skipped_lane = []
    kept = []

    for fq_path in fastqs:
        # Get basename from full GCS path
        basename = fq_path.split('/')[-1]

        if basename.startswith('Undetermined_') and not include_undetermined:
            skipped_undetermined.append(fq_path)
            continue

        sample_basename = None
        lane_num = None
        read_num = None

        match = pattern_lane.search(basename)
        if match:
            sample_basename = match.group(1)
            lane_num = int(match.group(2))
            read_num = match.group(3)
        else:
            match = pattern_nolane.search(basename)
            if match:
                sample_basename = match.group(1)
                read_num = match.group(2)

        # Only a *known* lane can disagree with the filter. A lane-less name means the
        # lanes were merged (or there is only one), so there is nothing to filter on.
        if lane_filter is not None and lane_num is not None and lane_num != lane_filter:
            skipped_lane.append(fq_path)
            continue

        # Inclusion is deliberately decoupled from parseability: anything surviving the
        # filters above lands in `fastqs`, even if its name did not parse. The consumer
        # downstream (tasks_demux.group_fastq_pairs) is more forgiving than this parser,
        # so dropping a file here would lose one it could have handled.
        kept.append(fq_path)

        if sample_basename is None:
            unparsed.append(fq_path)
            continue

        # Key on lane as well as sample, so a multi-lane run does not collapse its
        # per-lane files onto one another.
        group_key = sample_basename if lane_num is None else '{}_L{:03d}'.format(sample_basename, lane_num)
        if group_key not in groups:
            groups[group_key] = {}

        groups[group_key][f'R{read_num}'] = fq_path

    # Write kept fastqs to output
    with open('all_fastqs_output.txt', 'wt') as f:
        for fq in kept:
            f.write(fq + '\n')

    # Create output array of arrays
    # Each inner array is either [R1, R2] for PE or [R1] for SE
    pairs = []

    for group_key in sorted(groups.keys()):
        sample_files = groups[group_key]

        if 'R1' in sample_files and 'R2' in sample_files:
            # Paired-end
            pairs.append([sample_files['R1'], sample_files['R2']])
        elif 'R1' in sample_files:
            # Single-end (R1 only)
            pairs.append([sample_files['R1']])
        elif 'R2' in sample_files:
            # Single-end (R2 only, unusual but possible)
            pairs.append([sample_files['R2']])

    # Write output as JSON
    with open('raw_reads_fastq_pairs.json', 'wt') as f:
        json.dump(pairs, f, indent=2)

    def report(label, items):
        if not items:
            return
        print(f"{label}: {len(items)}")
        for item in items[:10]:  # Show first 10
            print(f"  {item}")
        if len(items) > 10:
            print(f"  ... and {len(items) - 10} more")

    report("Excluded Undetermined_* (set include_undetermined=true to keep)", skipped_undetermined)
    report(f"Excluded by lane filter (lane != {lane_filter})", skipped_lane)
    report("WARNING: included, but filename did not parse as *_R[12]_###.fastq.gz "
           "and so is absent from raw_reads_fastq_pairs", unparsed)

    # Fail loudly, and do it AFTER every filter so that an over-restrictive filter is
    # as loud as an empty directory. Returning an empty array under a green checkmark
    # is what kept this failure mode invisible.
    if not kept:
        print("ERROR: no usable FASTQ files.")
        print(f"  directory searched:        {fastq_dir_searched}")
        print(f"  .fastq.gz discovered:      {len(fastqs)}")
        print(f"  lane filter:               {lane_filter if lane_filter is not None else '(none)'}")
        print(f"  include_undetermined:      {include_undetermined}")
        print(f"  excluded as Undetermined:  {len(skipped_undetermined)}")
        print(f"  excluded by lane filter:   {len(skipped_lane)}")
        sys.exit(1)

    if lane_filter is not None:
        print(f"Filtered to lane {lane_filter}: {len(kept)} fastq files")
    print(f"Found {len(kept)} fastq files in {len(groups)} sample/lane groups with {len(pairs)} read groups")
    CODE
  >>>
  output {
    String               runinfo_xml           = read_string("runinfo_path.txt")
    Array[String]        fastqs                = read_lines("all_fastqs_output.txt")
    Array[Array[String]] raw_reads_fastq_pairs = read_json("raw_reads_fastq_pairs.json")
  }
  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: 1
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}
