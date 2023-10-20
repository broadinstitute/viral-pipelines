version 1.0

task gcs_copy {
  input {
    Array[File] infiles
    String      gcs_uri_prefix
  }
  meta {
    description: "gsutil cp. only works on a GCP-based backend (e.g. Terra)"
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
    gsutil -m cp ~{sep=' ' infiles} ~{gcs_uri_prefix}
  }
  output {
    File logs = stdout()
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-baseimage:0.1.20"
    memory: "1 GB"
    cpu: 1
    maxRetries: 1
  }
}

task get_gcloud_env_info {
  input {
    String? workspace_name
    String? workspace_namespace
    String? workspace_googleProject
  }
  meta {
    description: "task for inspection of backend env, and optionally running a command"
  }
  command <<<
    set -ex

    # ========== method 1: from google project IAM label metadata
    gcloud auth print-access-token

    GOOGLE_PROJECT_ID="$(gcloud config list --format='value(core.project)')"
    echo "$GOOGLE_PROJECT_ID" > google_project_id.txt

    #get project namespace:
    WORKSPACE_NAMESPACE="$(gcloud projects describe $GOOGLE_PROJECT_ID --format='value(labels.workspacenamespace)')"
    WORKSPACE_NAME="$(gcloud projects describe $GOOGLE_PROJECT_ID --format='value(labels.workspacename)')"

    echo "GOOGLE_PROJECT_ID:   ${GOOGLE_PROJECT_ID}"
    echo "WORKSPACE_NAMESPACE: ${WORKSPACE_NAMESPACE}"
    echo "WORKSPACE_NAME:      ${WORKSPACE_NAME}"

    # ========== method 2: matching a project returned by the API based on the google project ID

    GOOGLE_PROJECT_ID="$(gcloud config list --format='value(core.project)')"

    # get list of workspaces, limiting the output to only the fields we need
    curl -s -X 'GET' \
    'https://api.firecloud.org/api/workspaces?fields=workspace.name%2Cworkspace.namespace%2Cworkspace.googleProject' \
    -H 'accept: application/json' \
    -H "Authorization: Bearer $(gcloud auth print-access-token)" > workspace_list.json

    # extract workspace name
    WORKSPACE_NAME=$(jq -cr '.[] | select( .workspace.googleProject == "'${GOOGLE_PROJECT_ID}'" ).workspace | .name' workspace_list.json)
    echo "$WORKSPACE_NAME" > workspace_name.txt
    
    # extract workspace namespace
    WORKSPACE_NAMESPACE=$(jq -cr '.[] | select( .workspace.googleProject == "'${GOOGLE_PROJECT_ID}'" ).workspace | .namespace' workspace_list.json)
    echo "$WORKSPACE_NAMESPACE" > workspace_namespace.txt

    # ========== method 3: resolved by Terra as inputs

    echo ""
    echo "Strings passed in to workflow:"
    echo "workspace.name: ~{workspace_name}"
    echo "workspace.namespace: ~{workspace_namespace}"
    echo "workspace.googleProject: ~{workspace_googleProject}"
    echo ""

    # ======================================================================================

    env | tee -a env_info.log
    
    gcloud config list | tee -a gcloud_config_info.log
    
    gcloud info | tee -a gcloud_env_info.log

  >>>
  output {
    Array[String] env_info_files   = glob("./*_info.log")
    
    String workspace_name_out      = read_string("workspace_name.txt")
    String workspace_namespace_out = read_string("workspace_namespace.txt")
    String google_project_id_out   = read_string("google_project_id.txt")
  }
  runtime {
    docker: "quay.io/broadinstitute/viral-core:2.2.2"
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
    with open(out_fname, 'wt') as outf:
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
