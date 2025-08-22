version 1.0

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
  command {
    set -e
    gcloud storage cp "~{sep='" "' infiles}" ~{gcs_uri_prefix}
  }
  output {
    File logs = stdout()
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-baseimage:0.2.4"
    memory: "1 GB"
    cpu: 1
    maxRetries: 1
  }
}

task check_terra_env {
  input {
    String docker = "quay.io/broadinstitute/viral-core:2.4.1-23-g243eaaff-ct-swiftseq-demux-integration" #skip-global-version-pin
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
    maxRetries: 1
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
  command {
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
  }
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

    String        docker = "schaluvadi/pathogen-genomic-surveillance:api-wdl"
  }
  meta {
    volatile: true
  }
  command {
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
  }
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

    String  docker = "schaluvadi/pathogen-genomic-surveillance:api-wdl"
  }

  meta {
    volatile: true
  }

  command {
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
  }
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

    String  docker = "quay.io/broadinstitute/viral-core:2.4.1"
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
    lib_col_name      = "entity:~{library_table_name}_id"

    # import required packages
    import sys
    import collections
    import json
    import csv
    import pandas as pd
    import numpy as np
    from firecloud import api as fapi

    print(workspace_project + "\n" + workspace_name)

    # process read counts if available
    read_counts_raw = {}
    read_counts_cleaned = {}
    if '~{default="" read_counts_raw_json}':
        with open('~{default="" read_counts_raw_json}','rt') as inf:
            read_counts_raw = {pair['left']: pair['right'] for pair in json.load(inf)}
    if '~{default="" read_counts_cleaned_json}':
        with open('~{default="" read_counts_cleaned_json}','rt') as inf:
            read_counts_cleaned = {pair['left']: pair['right'] for pair in json.load(inf)}

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
    header, rows = get_entities_to_table(workspace_project, workspace_name, "~{sample_table_name}")
    df_sample = pd.DataFrame.from_records(rows, columns=header, index="~{sample_table_name}_id")
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
        outf.write('entity:~{sample_table_name}_id\tlibraries\n')
        merged_sample_ids = set()
        for sample_id, libraries in sample_to_libraries.items():
            if sample_id in df_sample.index and "libraries" in df_sample.columns and test_non_empty_value(df_sample.libraries[sample_id]):
                # merge in new sample->library mappings with any pre-existing sample->library mappings
                already_associated_libraries = [entity["entityName"] for entity in df_sample.libraries[sample_id] if entity.get("entityName")]
                libraries = list(set(libraries + already_associated_libraries))
                print (f"\tsample {sample_id} pre-exists in Terra table, merging old members {already_associated_libraries} with new members {libraries}")
                merged_sample_ids.add(sample_id)

            outf.write(f'{sample_id}\t{json.dumps([{"entityType":"~{library_table_name}","entityName":library_name} for library_name in libraries])}\n')
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
