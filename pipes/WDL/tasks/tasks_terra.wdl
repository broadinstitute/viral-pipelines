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

task check_terra_env {
  input {
    String docker = "quay.io/broadinstitute/viral-core:2.2.2" #skip-global-version-pin
  }
  meta {
    description: "task for inspection of backend to determine whether the task is running on Terra and/or GCP"
  }
  command <<<
    set -ex

    # create gcloud-related output file
    touch gcloud_config_info.log
    touch google_project_id.txt

    # create Terra-related output files
    touch workspace_name.txt
    touch workspace_namespace.txt
    touch workspace_bucket_path.txt
    touch input_table_name.txt
    touch input_row_id.txt

    # write system environment variables to output file
    env | tee -a env_info.log

    GOOGLE_PROJECT_ID="$(gcloud config list --format='value(core.project)')"
    echo "$GOOGLE_PROJECT_ID" > google_project_id.txt

    # check whether gcloud project has a "terra-" prefix
    # to determine if running on Terra
    if case ${GOOGLE_PROJECT_ID} in terra-*) ;; *) false;; esac; then
      # (shell-portable regex conditional)
      echo "Job appears to be running on Terra (GCP project ID: ${GOOGLE_PROJECT_ID})"
      echo "true" > RUNNING_ON_TERRA
    else
      echo "NOT running on Terra"
      echo "false" > RUNNING_ON_TERRA
    fi

    # check if running on GCP
    if curl -s metadata.google.internal -i | grep -E 'Metadata-Flavor:\s+Google'; then 
      echo "Cloud platform appears to be GCP"; 
      echo "true" > RUNNING_ON_GCP

      # write gcloud env info to output files
      gcloud info | tee -a gcloud_config_info.log
    else 
      echo "NOT running on GCP";
      echo "false" > RUNNING_ON_GCP
    fi

    if grep "true" RUNNING_ON_GCP && grep "true" RUNNING_ON_TERRA; then 
      echo "Running on Terra+GCP"

      # === Determine Terra workspace ID and submission ID for the workspace responsible for this job

      # Scrape various workflow / workspace info from the localization and delocalization scripts.
      #   from: https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/wdl/GvsUtils.wdl#L35-L40
      WORKSPACE_ID="$(sed -n -E 's!.*gs://fc-(secure-)?([^\/]+).*!\2!p' /cromwell_root/gcs_delocalization.sh | sort -u | tee workspace_id.txt)"
      echo "WORKSPACE_ID: ${WORKSPACE_ID}"

      # bucket path prefix
      #BUCKET_PREFIX="$(sed -n -E 's!.*(gs://(fc-(secure-)?[^\/]+)).*!\1!p' /cromwell_root/gcs_delocalization.sh | sort -u | tee bucket_prefix.txt)"
      #echo "BUCKET_PREFIX: ${BUCKET_PREFIX}"

      # top-level submission ID
      TOP_LEVEL_SUBMISSION_ID="$(sed -n -E 's!.*gs://fc-(secure-)?([^\/]+)/submissions/([^\/]+).*!\3!p' /cromwell_root/gcs_delocalization.sh | sort -u | tee top_level_submission_id.txt)"
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
        -H "Authorization: Bearer $(gcloud auth print-access-token)" > workspace_info.json


      WORKSPACE_NAME="$(jq -cr '.workspace.name | select (.!=null)' workspace_info.json)"
      WORKSPACE_NAME_URL_ENCODED="$(jq -rn --arg x "${WORKSPACE_NAME}" '$x|@uri')"
      WORKSPACE_NAMESPACE="$(jq -cr '.workspace.namespace | select (.!=null)' workspace_info.json)"
      WORKSPACE_BUCKET="gs://${WORKSPACE_ID}"

      echo "${WORKSPACE_NAME}" | tee workspace_name.txt
      echo "${WORKSPACE_NAMESPACE}" | tee workspace_namespace.txt
      echo "${WORKSPACE_BUCKET}" | tee workspace_bucket_path.txt

          # --- less direct way of obtaining workspace info by matching Terra project ID --
          #     preserved here for potential utility in obtaining workspace info for other projects/workspaces
          # get list of workspaces, limiting the output to only the fields we need
          #curl -s -X 'GET' \
          #'https://api.firecloud.org/api/workspaces?fields=workspace.name%2Cworkspace.namespace%2Cworkspace.bucketName%2Cworkspace.googleProject' \
          #-H 'accept: application/json' \
          #-H "Authorization: Bearer $(gcloud auth print-access-token)" > workspace_list.json

          # extract workspace name
          #WORKSPACE_NAME=$(jq -cr '.[] | select( .workspace.googleProject == "'${GOOGLE_PROJECT_ID}'" ).workspace | .name' workspace_list.json)
          
          # extract workspace namespace
          #WORKSPACE_NAMESPACE=$(jq -cr '.[] | select( .workspace.googleProject == "'${GOOGLE_PROJECT_ID}'" ).workspace | .namespace' workspace_list.json)
          #WORKSPACE_NAME_URL_ENCODED="$(jq -rn --arg x "${WORKSPACE_NAME}" '$x|@uri')"

          # extract workspace bucket
          #WORKSPACE_BUCKET=$(jq -cr '.[] | select( .workspace.googleProject == "'${GOOGLE_PROJECT_ID}'" ).workspace | .bucketName' workspace_list.json)
          # --- end less direct way of obtaining workspace info ---
      # =======================================


      # === obtain info on job submission inputs (table name, row ID)===
      touch submission_metadata.json
      curl -s -X 'GET' \
      "https://api.firecloud.org/api/workspaces/${WORKSPACE_NAMESPACE}/${WORKSPACE_NAME_URL_ENCODED}/submissions/${TOP_LEVEL_SUBMISSION_ID}" \
      -H 'accept: application/json' \
      -H "Authorization: Bearer $(gcloud auth print-access-token)" > submission_metadata.json

      INPUT_TABLE_NAME="$(jq -cr 'if .submissionEntity == null then "" elif (.workflows | length)==1 then .submissionEntity.entityType else [.workflows[].workflowEntity.entityType] | join(",") end' submission_metadata.json)"
      INPUT_ROW_ID="$(jq -cr 'if .submissionEntity == null then "" elif (.workflows | length)==1 then .submissionEntity.entityName else [.workflows[].workflowEntity.entityName] | join(",") end' submission_metadata.json)"

      echo "$INPUT_TABLE_NAME" | tee input_table_name.txt
      echo "$INPUT_ROW_ID" | tee input_row_id.txt
      # =======================================
    else 
      echo "Not running on Terra+GCP"
    fi

  >>>
  output {
    Boolean is_running_on_terra    = read_boolean("RUNNING_ON_TERRA")
    Boolean is_backed_by_gcp       = read_boolean("RUNNING_ON_GCP")

    String google_project_id       = read_string("google_project_id.txt")

    String workspace_id            = read_string("workspace_id.txt")
    String workspace_name          = read_string("workspace_name.txt")
    String workspace_namespace     = read_string("workspace_namespace.txt")
    String workspace_bucket_path   = read_string("workspace_bucket_path.txt")    

    String input_table_name        = read_string("input_table_name.txt")
    String input_row_id            = read_string("input_row_id.txt")

    String top_level_submission_id = read_string("top_level_submission_id.txt")

    File env_info                  = "env_info.log"
    File gcloud_config_info        = "gcloud_config_info.log"
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

task create_or_update_sample_tables {
  input {
    String flowcell_run_id

    String workspace_namespace
    String workspace_name
    String workspace_bucket

    String  docker = "quay.io/broadinstitute/viral-core:2.2.4" #skip-global-version-pin
  }

  meta {
    volatile: true
  }

  command <<<
    python3<<CODE

    flowcell_data_id   = '~{flowcell_run_id}'
    
    workspace_project = '~{workspace_namespace}'
    workspace_name    = '~{workspace_name}'
    workspace_bucket  = '~{workspace_bucket}'

    # import required packages.
    import os
    import argparse
    import collections
    import json
    import csv
    import pandas as pd
    from firecloud import api as fapi
    from ast import literal_eval
    from io import StringIO

    print(workspace_project + "\n" + workspace_name + "\n" + "bucket: " + workspace_bucket)

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

    # # populate sample table from run outputs
    # In particular, popualte the raw_bam and cleaned_bam columns

    def create_tsvs(project, workspace, runID):
        # API call to get flowcell_data table
        response = fapi.get_entities_tsv(project, workspace, "flowcell", model="flexible")

        # read API response into data frame
        df = pd.read_csv(StringIO(response.text), sep="\t", index_col="entity:flowcell_id")

        # create sample.tsv data frame (entity:sample_set_id)
        cleaned_bams_list = literal_eval(df.cleaned_reads_unaligned_bams[runID])
        cleaned_library_id_list = [bam.split("/")[-1].replace(".bam", "").replace(".cleaned", "") for bam in cleaned_bams_list]
        df_library_table = pd.DataFrame({"entity:library_id" : cleaned_library_id_list,
                                        "cleaned_bam" : cleaned_bams_list})
        cleaned_library_fname = runID + "_cleaned_bams.tsv"
        df_library_table.to_csv(cleaned_library_fname, sep="\t", index=False)

        # create sample.tsv data frame (entity:sample_set_id)
        raw_bams_list = literal_eval(df.raw_reads_unaligned_bams[runID])
        raw_library_id_list = [bam.split("/")[-1].replace(".bam", "") for bam in raw_bams_list]
        df_library_table = pd.DataFrame({"entity:library_id" : raw_library_id_list,
                                        "raw_bam" : raw_bams_list})
        raw_library_fname = runID + "_raw_bams.tsv"
        df_library_table.to_csv(raw_library_fname, sep="\t", index=False)

        ## create sample_set.tsv data frame (membership:sample_set_id)
        #participant_id_list = [sample_id.split(".")[0] for sample_id in raw_bams_list]
        #df_sample_set_table = pd.DataFrame({"membership:sample_set_id" : participant_id_list,
        #                                    "sample" : sample_id_list})
        #set_fname = runID + "_sample_set_table.tsv"
        #df_sample_set_table.to_csv(set_fname, sep="\t", index=False)
        
        return (cleaned_library_fname, raw_library_fname)

    # call the create_tsv function and save files to data tables
    tsv_list = create_tsvs(workspace_project, workspace_name, flowcell_data_id)

    print (f"wrote outputs to {tsv_list}")

    for tsv in tsv_list:
        response = fapi.upload_entities_tsv(workspace_project, workspace_name, tsv, model="flexible")
        if response.status_code != 200:
            print('ERROR UPLOADING: See full error message:')
            print(response.text)
        else:
            print("Upload complete. Check your workspace for new table!")

    # # update sample_set with new set memberships and flowcell metadata

    # columns to copy from flowcell_data to library table
    copy_cols = ["sample_original", "spike_in"]

    # API call to get existing sample_set mappings
    header, rows = get_entities_to_table(workspace_project, workspace_name, "sample")
    df_sample = pd.DataFrame.from_records(rows, columns=header, index="sample_id")

    # API call to get all existing library ids
    header, rows = get_entities_to_table(workspace_project, workspace_name, "library")
    df_library = pd.DataFrame.from_records(rows, columns=header, index="library_id")

    # API call to get flowcell_data table
    header, rows = get_entities_to_table(workspace_project, workspace_name, "flowcell")
    df_flowcell = pd.DataFrame.from_records(rows, columns=header, index="flowcell_id")

    # grab the meta_by_filename values to create new sample->library (sample_set->sample) mappings
    sample_to_libraries = {}
    for library_id, data in df_flowcell.meta_by_filename[flowcell_data_id].items():
        sample_id = data['sample']
        sample_to_libraries.setdefault(sample_id, [])
        if library_id in df_library.index:
            sample_to_libraries[sample_id].append(library_id)
        else:
            print (f"missing {library_id} from sample table")

    # merge in new sample->library mappings with any pre-existing sample->library mappings
    if len(df_sample)>0:
        print(df_sample.index)
        for sample_id in sample_to_libraries.keys():
            if sample_id in df_sample.index:
                print (f"sample_set {sample_id} pre-exists in Terra table, merging with new members")
                #sample_set_to_samples[set_id].extend(df_sample_set.samples[set_id])
                already_associated_libraries = [entity["entityName"] for entity in df_sample.libraries[sample_id]]
                
                print(f"already_associated_libraries {already_associated_libraries}")
                print(f"sample_to_libraries[sample_id] {sample_to_libraries[sample_id]}")
                continue
                sample_to_libraries[sample_id].extend(already_associated_libraries)
                # collapse duplicate sample IDs
                #sample_to_libraries[sample_id] = list([json.dumps({"entityType":"library","entityName":library_name}) for library_name in set(sample_to_libraries[sample_id])])
                sample_to_libraries[sample_id] = list(set(sample_to_libraries[sample_id]))
                #sample_to_libraries[sample_id] = list({"entityType":"library","entityName":library_name}set(sample_to_libraries[sample_id]))

    # create sample_membership.tsv describing sample:library (ie sample-to-library) many-to-one mappings
    # sample_fname = 'sample_membership.tsv'
    # with open(sample_fname, 'wt') as outf:
    #     outf.write('membership:sample_id\tlibrary\n')
    #     for sample_id, libraries in sample_to_libraries.items():
    #         for library_id in sorted(libraries):
    #             outf.write(f'{sample_id}\t{library_id}\n')
    # !cat $sample_fname

    sample_fname = 'sample_membership.tsv'
    with open(sample_fname, 'wt') as outf:
        outf.write('entity:sample_id\tlibraries\n')
        for sample_id, libraries in sample_to_libraries.items():
            #for library_id in sorted(libraries):
            outf.write(f'{sample_id}\t{json.dumps([{"entityType":"library","entityName":library_name} for library_name in libraries])}\n')

    # grab the meta_by_sample values from one row in the flowcell_data table
    meta_by_library_all = df_flowcell.meta_by_sample[flowcell_data_id]

    # grab all the library IDs
    header, rows = get_entities_to_table(workspace_project, workspace_name, "library")
    out_rows = []
    out_header = ['library_id'] + copy_cols
    print(f"out_header {out_header}")
    for row in rows:
        out_row = {'library_id': row['library_id']}

        for sample_id,sample_library_metadata in meta_by_library_all.items():
            if sample_library_metadata["library"] in row['library_id']:
                for col in copy_cols:
                    out_row[col] = sample_library_metadata.get(col, '')
                out_rows.append(out_row)

    library_meta_fname = "sample_metadata.tsv"
    with open(library_meta_fname, 'wt') as outf:
      outf.write("entity:")
      writer = csv.DictWriter(outf, out_header, delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
      writer.writeheader()
      writer.writerows(out_rows)

    # write them to the Terra table!
    for fname in (library_meta_fname,sample_fname):
        response = fapi.upload_entities_tsv(workspace_project, workspace_name, fname, model="flexible")
        if response.status_code != 200:
            print(f'ERROR UPLOADING {fname}: See full error message:')
            print(response.text)
        else:
            print("Upload complete. Check your workspace for new table!")

    CODE
  >>>
  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    maxRetries: 2
  }
  output {
    File stdout_log = stdout()
    File stderr_log = stderr()
  }
}
