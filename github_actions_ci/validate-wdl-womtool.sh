#!/bin/bash
set -e -o pipefail

# if WORKFLOW_TO_VALIDATE env var is set, validate only that workflow
#   ex. WORKFLOW_TO_VALIDATE=demux_deplete ./github_actions_ci/validate-wdl-womtool.sh

if [ -n "$WORKFLOW_TO_VALIDATE" ]; then
  WF_GLOB_PATTERN="pipes/WDL/workflows/${WORKFLOW_TO_VALIDATE}.wdl"
else
  WF_GLOB_PATTERN='pipes/WDL/workflows/*.wdl'
fi

for workflow in ${WF_GLOB_PATTERN}; do
  printf "validating $workflow ... "
  if $(hash -r  womtool &> /dev/null); then
    womtool validate $workflow
  else
    java -jar womtool.jar validate $workflow
  fi
done
