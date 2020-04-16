#!/bin/bash
set -ex -o pipefail

mkdir -p workflows
cp -r test workflows/
cd workflows

miniwdl run_self_test

for workflow in ../pipes/WDL/workflows/*.wdl; do
	workflow_name=$(basename $workflow .wdl)
	input_json="test/input/WDL/test_inputs-$workflow_name-local.json"
	if [ -f $input_json ]; then
		date
		echo "validating $workflow"
		miniwdl check $workflow
		echo "Executing $workflow_name using miniWDL on local instance"
		miniwdl run -i $input_json -d $wfname/. $workflow
		if [ -f $workflow_name/outputs.json ]; then
			echo "$workflow_name SUCCESS -- outputs:"
			cat $workflow_name/outputs.json
		else
			echo "$workflow_name FAILED"
			exit 1
		fi
    fi
done

cd -
date
echo "note: there is no testing of output correctness yet..."
