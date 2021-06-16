#!/bin/bash
set -e  # intentionally allow for pipe failures below

mkdir -p workflows
cp *.jar pipes/WDL/workflows/*.wdl pipes/WDL/tasks/*.wdl workflows
cp -r test workflows/
cd workflows

for workflow in ../pipes/WDL/workflows/*.wdl; do
	workflow_name=$(basename $workflow .wdl)
	input_json="test/input/WDL/cromwell-local/test_inputs-$workflow_name-local.json"
	if [ -f $input_json ]; then
		date
		echo "Executing $workflow_name using Cromwell on local instance"
		# the "cat" is to allow a pipe failure (otherwise it halts because of set -e)
		java -Dconfig.file=../pipes/cromwell/cromwell.local-github_actions.conf \
			-jar cromwell.jar run \
			$workflow_name.wdl \
			-i $input_json | tee cromwell.out
		if [ ${PIPESTATUS[0]} -gt 0 ]; then
			echo "error running $workflow_name"
			error_logs=$(grep stderr cromwell.out | perl -lape 's/.*\s(\S+)$/$1/g')
			for log in $error_logs; do
				echo "contents of stderr ($log):"
				cat `dirname $log`/stderr | sed "s/^/[STDERR] /"
				echo "contents of stdout ($log):"
				cat `dirname $log`/stdout | sed "s/^/[STDOUT] /"
			done
			sync; sleep 20; exit 1
		fi
    fi
done

cd -
date
echo "note: there is no testing of output correctness yet..."
