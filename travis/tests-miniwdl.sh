#!/bin/bash
set -e -o pipefail  # intentionally allow for pipe failures below

mkdir -p workflows
ln -s test workflows/
cd workflows

miniwdl run_self_test

for workflow in ../pipes/WDL/workflows/*.wdl; do
	workflow_name=$(basename $workflow .wdl)
	input_json="test/input/WDL/test_inputs-$workflow_name-local.json"
	if [ -f $input_json ]; then
		date
		echo "validating $workflow"
		wfname=`basename $workflow`
		miniwdl check $workflow
		echo "Executing $workflow_name using miniWDL on local instance"
		miniwdl run -i $input_json -d $wfname/. $workflow
		if [ -f $wfname/outputs.json ]; then
			echo "$wfname SUCCESS -- outputs:"
			cat $wfname/outputs.json
		else
			echo "$wfname FAILED"
			exit 1
		fi
    fi
done

cd -
date
echo "note: there is no testing of output correctness yet..."
