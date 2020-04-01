#!/bin/bash
set -e -o pipefail

FLAT_DIR=../pipes/WDL/flattened
mkdir -p $FLAT_DIR

for workflow in ../pipes/WDL/workflows/*.wdl; do
	wf_base=`basename $workflow`
	out_fn="$FLAT_DIR/$wf_base"
	echo "flattening $workflow to $out_fn"
	./paste_wdl_imports.py -o $out_fn $workflow

	### DEBUG
	cat $out_fn
done
