#!/bin/bash
set -e -o pipefail

for workflow in pipes/WDL/workflows/*.wdl; do
  echo -n "validating $workflow ... "
  if $(hash -r  womtool &> /dev/null); then
    womtool validate $workflow
  else
    java -jar womtool.jar validate $workflow
  fi
done
