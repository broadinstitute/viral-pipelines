#!/bin/bash
set -e -o pipefail

for workflow in ../pipes/WDL/workflows/*.wdl; do
  echo "validating $workflow"
  miniwdl check $workflow
done
