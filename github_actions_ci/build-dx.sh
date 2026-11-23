#!/bin/bash
set -e -o pipefail

if [ -z "$DX_API_TOKEN" ]; then
  echo "ERROR: DX_API_TOKEN is not set, this is needed to build DNAnexus workflows."
  exit 1
fi

# obtain version tag
VERSION=`github_actions_ci/list-docker-tags.sh | tail -1 | sed 's/:/\//'`

# log in to DNAnexus
source dx-toolkit/environment
dx login --token "$DX_API_TOKEN" --noprojects
dx select $DX_PROJECT

# compile with dxCompiler
COMPILE_SUCCESS="dxCompiler-compile_all-success.txt"
DX_COMPILER_PARALLELISM="${DX_COMPILER_PARALLELISM:-4}"
DX_COMPILER_TIMEOUT="${DX_COMPILER_TIMEOUT:-15m}"
echo "dxCompiler executable creation parallelism: $DX_COMPILER_PARALLELISM"
echo "dxCompiler per-workflow timeout: $DX_COMPILER_TIMEOUT"
touch $COMPILE_SUCCESS
for workflow in pipes/WDL/workflows/*.wdl; do
  if [ -n "$(grep DX_SKIP_WORKFLOW $workflow)" ]; then
    echo "Skipping $workflow due to the presence of the DX_SKIP_WORKFLOW tag"
  else
    workflow_name=`basename $workflow .wdl`
	  echo "Building $workflow to DNAnexus: /build/$VERSION/$workflow_name"

    defaults_json="pipes/dnax/dx-defaults-$workflow_name.json"
    if [ -f "$defaults_json" ]; then
      CMD_DEFAULTS=(-defaults "$defaults_json")
    else
      CMD_DEFAULTS=()
    fi

    extras_json="pipes/dnax/dx-extras.json"
    CMD_DEFAULTS+=(-extras "$extras_json")

    compile_start=$SECONDS
    if dx_id=$(timeout --signal=TERM --kill-after=30s "$DX_COMPILER_TIMEOUT" \
      java -jar dxCompiler.jar compile \
        "$workflow" "${CMD_DEFAULTS[@]}" -f \
        -executableCreationParallelism "$DX_COMPILER_PARALLELISM" \
        -leaveWorkflowsOpen \
        -imports pipes/WDL/tasks/ \
        -project "$DX_PROJECT" \
        -destination "/build/$VERSION/$workflow_name"); then
        echo "Succeeded in $((SECONDS - compile_start))s: $workflow_name = $dx_id"
    else
        compile_status=$?
        if [ "$compile_status" -eq 124 ]; then
          echo "Timed out after $((SECONDS - compile_start))s building: $workflow_name"
        else
          echo "Failed after $((SECONDS - compile_start))s building: $workflow_name (exit $compile_status)"
        fi
        exit "$compile_status"
    fi
    echo -e "$workflow_name\t$dx_id" >> $COMPILE_SUCCESS
  fi
done

# Special case: build demux launchers (native DNAnexus applets), embedding the
# demux workflow ID as a default input. Skip if no demux workflows were compiled.
demux_workflows_to_build="demux_plus demux_only"
any_demux_compiled=false
for wf_name in $(echo "${demux_workflows_to_build}"); do
  if grep -q "^${wf_name}\s" $COMPILE_SUCCESS; then
    any_demux_compiled=true
  fi
done

if [ "$any_demux_compiled" = true ]; then
  # build consolidate_run_tarballs (native DNAnexus applet) applet
  pushd pipes/dnax/dx-launcher
  cp consolidate_run_tarballs.yml consolidate_run_tarballs_dxapp.yml
  consolidate_tarballs_dx_id=$(./dx-yml-build consolidate_run_tarballs_dxapp.yml -a --destination /build/$VERSION/ | jq -r ".id")
  popd
  echo -e "consolidate_run_tarballs\t$consolidate_tarballs_dx_id" >> $COMPILE_SUCCESS

  for wf_name in $(echo "${demux_workflows_to_build}"); do
    demux_workflow_id=$(grep "^${wf_name}\s" $COMPILE_SUCCESS | cut -f 2)
    if [ -z "$demux_workflow_id" ]; then
      echo "Skipping applet ${wf_name}_launcher (${wf_name} was not compiled)"
      continue
    fi
    echo "Building applet ${wf_name}..."
    pushd pipes/dnax/dx-launcher
    sed "s/DEFAULT_DEMUX_WORKFLOW_ID/$demux_workflow_id/" demux_launcher.yml \
      | sed "s/DEFAULT_DEMUX_WORKFLOW_NAME/${wf_name}_launcher/" \
      | sed "s/DEFAULT_CONSOLIDATE_RUN_TARBALLS_APPLET_ID/$consolidate_tarballs_dx_id/" > "${wf_name}_dxapp.yml"
    dx_id=$(./dx-yml-build ${wf_name}_dxapp.yml -a --destination /build/$VERSION/ | jq -r ".id")
    popd
    echo -e "${wf_name}_launcher\t$dx_id" >> $COMPILE_SUCCESS
  done
else
  echo "Skipping consolidate_run_tarballs and demux launchers (no demux workflows were compiled)"
fi

# the presence of this file in the project denotes successful build
dx upload --brief --no-progress --destination /build/$VERSION/ $COMPILE_SUCCESS
