#!/bin/bash
set -ex -o pipefail

starting_dir="$(pwd)"
test_dir="miniwdl_testing"

function cleanup(){
    echo "Cleaning up from miniwdl run; exit code: $?"
    cd "$starting_dir"
    if [ -d "$test_dir" ]; then
      rm -r "$test_dir"
    fi
}
trap cleanup EXIT SIGINT SIGQUIT SIGTERM

mkdir $test_dir
cp -r test $test_dir
cd $test_dir

docker --version

# make sure our system has everything it needs to perform "miniwdl run" (e.g. docker swarm works)
miniwdl run_self_test

for workflow in ../pipes/WDL/workflows/*.wdl; do
	workflow_name=$(basename $workflow .wdl)
	input_json="test/input/WDL/miniwdl-local/test_inputs-$workflow_name-local.json"
	expected_output_json="test/input/WDL/miniwdl-local/test_outputs-$workflow_name-local.json"
	if [ -f $input_json ]; then
		date
		echo "Executing $workflow_name using miniWDL on local instance"
		# the following invocation with -d $workflow_name/. tells miniwdl not to create
		# a timestamped subdirectory within $workflow_name/
		time miniwdl run -i $input_json -d $workflow_name/. --error-json --verbose $workflow
		# the existence of $workflow_name/outputs.json is a guarantee of successful execution
		if [ -f $workflow_name/outputs.json ]; then
			echo "$workflow_name SUCCESS -- outputs:"
			cat $workflow_name/outputs.json
			if [ -f $expected_output_json ]; then
				echo "$workflow_name -- validating outputs"
				touch expected actual
				# here we create a key=value text table for the output in "actual" for only the
				# subset of keys that are specified in $expected_output_json, ignoring any
				# of the actual output values that aren't specified there.
				for k in `cat $expected_output_json | jq -r 'keys[]'`; do
					echo -n "$k=" >> expected
					echo -n "$k=" >> actual
					cat $expected_output_json       | jq -r '.["'$k'"]' >> expected
					cat $workflow_name/outputs.json | jq -r '.["'$k'"]' >> actual
				done
				# diff returns non-zero (and triggers failure by set -e) if these are non-equal
				# and prints all differences between expected and actual for this workflow at the same time
				diff expected actual
			fi
		else
			echo "$workflow_name FAILED"
			exit 1
		fi
		docker image prune --all --force # prune images from this workflow execution to save space before the next
	fi
done

cd "$starting_dir"
date
