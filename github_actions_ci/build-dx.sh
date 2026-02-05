#!/bin/bash
set -e -o pipefail

if [ -z "$DX_API_TOKEN" ]; then
  echo "ERROR: DX_API_TOKEN is not set, this is needed to build DNAnexus workflows."
  exit 1
fi

# Keep-alive mechanism: print status every 60s, abort if no output for 20 mins
# This prevents GitHub Actions from killing the job due to no output
SILENCE_TIMEOUT=1200  # 20 minutes in seconds
LAST_OUTPUT_FILE=$(mktemp)
date +%s > "$LAST_OUTPUT_FILE"

# Background keep-alive process
(
  while true; do
    sleep 60
    LAST_OUTPUT=$(cat "$LAST_OUTPUT_FILE" 2>/dev/null || echo 0)
    NOW=$(date +%s)
    SILENT_FOR=$((NOW - LAST_OUTPUT))
    echo "$(date '+%Y-%m-%d %H:%M:%S'): dxCompiler running (${SILENT_FOR}s since last output)..."
    if [ "$SILENT_FOR" -gt "$SILENCE_TIMEOUT" ]; then
      echo "ERROR: No dxCompiler output for ${SILENCE_TIMEOUT}s, aborting"
      kill -TERM $$ 2>/dev/null || true
      exit 1
    fi
  done
) &
KEEPALIVE_PID=$!
trap "kill $KEEPALIVE_PID 2>/dev/null; rm -f $LAST_OUTPUT_FILE" EXIT

# Wrapper function to run commands while tracking output timestamps
run_with_keepalive() {
  # Run command, update timestamp on each line, suppress broken pipe errors
  # The 2>/dev/null on echo suppresses "Broken pipe" when the command substitution closes
  "$@" 2>&1 | while IFS= read -r line; do
    date +%s > "$LAST_OUTPUT_FILE"
    echo "$line" 2>/dev/null || true
  done
  # Return the exit status of the first command in the pipeline
  return ${PIPESTATUS[0]}
}

# obtain version tag
VERSION=`github_actions_ci/list-docker-tags.sh | tail -1 | sed 's/:/\//'`

# log in to DNAnexus
source dx-toolkit/environment
dx login --token "$DX_API_TOKEN" --noprojects
dx select $DX_PROJECT

# compile with dxCompiler
COMPILE_SUCCESS="dxCompiler-compile_all-success.txt"
touch $COMPILE_SUCCESS
for workflow in pipes/WDL/workflows/*.wdl; do
  if [ -n "$(grep DX_SKIP_WORKFLOW $workflow)" ]; then
    echo "Skipping $workflow due to the presence of the DX_SKIP_WORKFLOW tag"
  else
    workflow_name=`basename $workflow .wdl`
	  echo "Building $workflow to DNAnexus: /build/$VERSION/$workflow_name"

    defaults_json="pipes/dnax/dx-defaults-$workflow_name.json"
    if [ -f "$defaults_json" ]; then
      CMD_DEFAULTS="-defaults $defaults_json"
    else
      CMD_DEFAULTS=""
    fi

    extras_json="pipes/dnax/dx-extras.json"
    CMD_DEFAULTS+=" -extras $extras_json"

	  dx_id=$(run_with_keepalive java -jar dxCompiler.jar compile \
      $workflow $CMD_DEFAULTS -f -verbose \
      -leaveWorkflowsOpen \
      -imports pipes/WDL/tasks/ \
      -project $DX_PROJECT \
      -destination /build/$VERSION/$workflow_name)
    if [ $? -eq 0 ]; then
        echo "Succeeded: $workflow_name = $dx_id"
    else
        echo "Failed to build: $workflow_name"
        exit $?
    fi
    echo -e "$workflow_name\t$dx_id" >> $COMPILE_SUCCESS
  fi
done

# build consolidate_run_tarballs (native DNAnexus applet) applet
pushd pipes/dnax/dx-launcher
cp consolidate_run_tarballs.yml consolidate_run_tarballs_dxapp.yml
consolidate_tarballs_dx_id=$(./dx-yml-build consolidate_run_tarballs_dxapp.yml -a --destination /build/$VERSION/ | jq -r ".id")
popd
echo -e "consolidate_run_tarballs\t$dx_id" >> $COMPILE_SUCCESS

# Special case: build demux launchers (native DNAnexus applets), embedding the
# demux workflow ID as a default input 
demux_workflows_to_build="demux_plus demux_only"
for wf_name in $(echo "${demux_workflows_to_build}"); do
  echo "Building applet ${wf_name}..."
  demux_workflow_id=$(grep "^${wf_name}\s" $COMPILE_SUCCESS | cut -f 2)
  pushd pipes/dnax/dx-launcher
  sed "s/DEFAULT_DEMUX_WORKFLOW_ID/$demux_workflow_id/" demux_launcher.yml \
    | sed "s/DEFAULT_DEMUX_WORKFLOW_NAME/${wf_name}_launcher/" \
    | sed "s/DEFAULT_CONSOLIDATE_RUN_TARBALLS_APPLET_ID/$consolidate_tarballs_dx_id/" > "${wf_name}_dxapp.yml"
  dx_id=$(./dx-yml-build ${wf_name}_dxapp.yml -a --destination /build/$VERSION/ | jq -r ".id")
  popd
  echo -e "${wf_name}_launcher\t$dx_id" >> $COMPILE_SUCCESS
done

# the presence of this file in the project denotes successful build
dx upload --brief --no-progress --destination /build/$VERSION/ $COMPILE_SUCCESS
