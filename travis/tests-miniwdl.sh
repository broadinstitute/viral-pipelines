#!/bin/bash
set -ex -o pipefail

mkdir -p workflows
cp -r test workflows/
cd workflows

miniwdl run_self_test

for workflow in ../pipes/WDL/workflows/*.wdl; do
	workflow_name=$(basename $workflow .wdl)
	input_json="test/input/WDL/test_inputs-$workflow_name-local.json"
	expected_output_json="test/input/WDL/test_outputs-$workflow_name-local.json"
	if [ -f $input_json ]; then
		date
		echo "Executing $workflow_name using miniWDL on local instance"
		miniwdl run -i $input_json -d $workflow_name/. --error-json $workflow
		if [ -f $workflow_name/outputs.json ]; then
			echo "$workflow_name SUCCESS -- outputs:"
			cat $workflow_name/outputs.json
			if [ -f $expected_output_json ]; then
				echo "$workflow_name -- validating outputs"
				touch expected actual
				for k in `cat $expected_output_json | jq -r 'keys[]'`; do
					echo -n "$k=" >> expected
					echo -n "$k=" >> actual
					cat $expected_output_json       | jq -r '.["'$k'"]' >> expected
					cat $workflow_name/outputs.json | jq -r '.["'$k'"]' >> actual
				done
				diff expected actual
			fi
		else
			echo "$workflow_name FAILED"
			exit 1
		fi
	fi
done

cd -
date
