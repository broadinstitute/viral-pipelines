#!/bin/bash
set -e  # intentionally allow for pipe failures below

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
		miniwdl check $workflow
		echo "Executing $workflow_name using miniWDL on local instance"
		miniwdl run -i $input_json $workflow
    fi
done

cd -
date
echo "note: there is no testing of output correctness yet..."
