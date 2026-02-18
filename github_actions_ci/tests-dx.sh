#!/bin/bash
set -e -o pipefail

# obtain version tag
VERSION=`github_actions_ci/list-docker-tags.sh | tail -1 | sed 's/:/\//'`

# log in to DNAnexus
source dx-toolkit/environment
if [ -n "$DX_API_TOKEN" ]; then
  dx login --token "$DX_API_TOKEN" --noprojects
  dx select $DX_PROJECT
fi

COMPILE_SUCCESS="dxCompiler-compile_all-success.txt"
if [ ! -f $COMPILE_SUCCESS ]; then
  dx download --no-progress /build/$VERSION/$COMPILE_SUCCESS
fi

function dx_run_timeout_args {
    #
    # Construct command-line arguments for 'dx run' command
    # to set a timeout on the applets it runs
    #

    local dx_workflow_id="$1"
    local dx_extra_applet_id="$2"

    local dx_workflow_applet_ids=$(dx describe $dx_workflow_id | grep applet- | awk '{print $2;}')
    local dx_applet_ids="$dx_workflow_applet_ids $dx_extra_applet_id"
    local comma=""
    local timeout_args="{\"timeoutPolicyByExecutable\":{"
    for dx_applet_id in $dx_applet_ids
    do
        timeout_args="${timeout_args}${comma}\"$dx_applet_id\":{\"*\":{\"hours\":3}}"
        comma=","
    done
    timeout_args="$timeout_args}}"
    echo $timeout_args
}

TEST_LAUNCH_ALL="dxWDL-execute_all-launched.txt"
touch $TEST_LAUNCH_ALL
for workflow in pipes/WDL/workflows/*.wdl; do
  echo "testing $workflow..."
  if [ -n "$(grep DX_SKIP_WORKFLOW $workflow)" ]; then
    echo "Skipping $workflow due to the presence of the DX_SKIP_WORKFLOW tag"
  else
    workflow_name=`basename $workflow .wdl`
    input_json="test/input/WDL/test_inputs-$workflow_name-dnanexus.dx.json"
    if [ -f $input_json ]; then
       # launch simple test cases on DNAnexus CI project
       dx_workflow_id=$(grep -w "^$workflow_name" $COMPILE_SUCCESS | cut -f 2)
       timeout_args=$(dx_run_timeout_args $dx_workflow_id)
       echo "running test $workflow_name - $dx_workflow_id -y --brief -f $input_json --extra-args $timeout_args"
       dx_job_id=$(dx run \
           $dx_workflow_id -y --brief \
           -f $input_json \
           --name "$VERSION $workflow_name" \
           --destination /tests/$VERSION/$workflow_name \
           #--extra-args $timeout_args \
           )
       if [ $? -eq 0 ]; then
           echo "Launched $workflow_name as $dx_job_id"
       else
           echo "Failed to build: $workflow_name"
       fi
       echo "Launched $workflow_name as $dx_job_id"
       echo -e "$workflow_name\t$dx_workflow_id\t$dx_job_id" >> $TEST_LAUNCH_ALL
    fi
  fi
done

# only run demux_plus if this is on master or tagged branch
if [ "$GITHUB_ACTIONS_BRANCH" = "master" -o -n "$GITHUB_ACTIONS_TAG" ]; then
  demux_name="demux_plus"
else
  # otherwise just run the (faster) demux_only
  demux_name="demux_only"
fi

# Special case: run test for the demux_(plus|only)_launcher native applet (which invokes
# the demux_(plus|only) WDL workflow). Skip if the launcher was not built.
demux_launcher_id=$(grep "^${demux_name}_launcher\s" $COMPILE_SUCCESS | cut -f 2 || true)
demux_workflow_id=$(grep "^${demux_name}\s" $COMPILE_SUCCESS | cut -f 2 || true)

if [ -n "$demux_launcher_id" -a -n "$demux_workflow_id" ]; then
  timeout_args=$(dx_run_timeout_args $demux_workflow_id $demux_launcher_id)
  dx_job_id=$(dx run "${demux_launcher_id}" \
    -y --brief \
    -i upload_sentinel_record=record-Bv8qkgQ0jy198GK0QVz2PV8Y \
    -i demux_workflow_id=${demux_workflow_id} \
    --name "$VERSION ${demux_name}_launcher" \
    -i folder=/tests/$VERSION/${demux_name}_launcher \
    --extra-args $timeout_args \
    )
  echo "Launched ${demux_name}_launcher as $dx_job_id"
  echo -e "${demux_name}_launcher\t$demux_launcher_id\t$dx_job_id" >> $TEST_LAUNCH_ALL
else
  echo "Skipping ${demux_name}_launcher test (launcher was not built)"
fi

# the presence of this file in the project denotes all tests launched
dx upload --brief --no-progress --destination /build/$VERSION/ $TEST_LAUNCH_ALL

#
# Cleanup folders w/ files that are 30 days old (@yihchii)
#
THIRTY_DAYS_FILES="30d_old_files.json"
dx find data --class file --created-before=-30d --path=/tests/ --json > $THIRTY_DAYS_FILES
[ "$(cat $THIRTY_DAYS_FILES)" != "[]" ] && dx rm -r -f $(jq -r '(.[].describe.folder + "/")' $THIRTY_DAYS_FILES | sort -u)
rm -rf $THIRTY_DAYS_FILES

# Cleanup empty folders
#  - Loop over all versioned sub-folders and clean them up if they are empty
ALL_VERSION_FOLDERS="all_version_folders.txt"
CURRENT_FOLDER="/tests/$VERSION/"
dx ls --folders --full $(dirname $CURRENT_FOLDER) > $ALL_VERSION_FOLDERS
for folder in $(cat $ALL_VERSION_FOLDERS); do 
    content=$(dx find data --path $folder --brief)
    [ -z "$content" ] && dx rm -rf $folder
done
rm -rf $ALL_VERSION_FOLDERS