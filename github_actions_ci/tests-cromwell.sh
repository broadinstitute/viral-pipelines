#!/bin/bash
set -e  # intentionally allow for pipe failures below

# increase docker timeouts to allow for staging of larger images (seconds)
export DOCKER_CLIENT_TIMEOUT=240
export COMPOSE_HTTP_TIMEOUT=240

starting_dir="$(pwd)"
test_dir="cromwell_testing"

function cleanup(){
    echo "Cleaning up from miniwdl run; exit code: $?"
    cd "$starting_dir"
    if [ -d "$test_dir" ] && [[ $KEEP_OUTPUT != "true" ]]; then
      rm -r "$test_dir"
    fi
}
trap cleanup EXIT SIGINT SIGQUIT SIGTERM

mkdir -p ${test_dir}
cp pipes/WDL/workflows/*.wdl pipes/WDL/tasks/*.wdl $test_dir
sed -i -- 's|import \"../tasks/|import \"|g' ${test_dir}/*.wdl
cp -r test ${test_dir}/
cd ${test_dir}

CROMWELL_LOG_LEVEL="${CROMWELL_LOG_LEVEL:=WARN}"

# if "cromwell" exists on the PATH (no .jar file extension suffix)
# it means it was installed from bioconda
if hash cromwell &>/dev/null; then
	echo "conda cromwell present";
	# this is the bioconda java-launching script
	JAVA_ENTRYPOINT="cromwell"
else
	# otherwise if cromwell is not installed via conda, call java
	JAVA_ENTRYPOINT="java"
	cp *.jar ${test_dir}
	CROMWELL_JAR_ARG="-jar cromwell.jar"
fi

for workflow in ../pipes/WDL/workflows/*.wdl; do
	workflow_name=$(basename $workflow .wdl)
	input_json="test/input/WDL/cromwell-local/test_inputs-$workflow_name-local.json"
	if [ -f $input_json ]; then
		date
		echo "Executing $workflow_name using Cromwell on local instance"
		# the "cat" is to allow a pipe failure (otherwise it halts because of set -e)
		${JAVA_ENTRYPOINT} -Dconfig.file=../pipes/cromwell/cromwell.local-github_actions.conf \
			-DLOG_MODE=pretty \
			-DLOG_LEVEL=${CROMWELL_LOG_LEVEL} \
			${CROMWELL_JAR_ARG} run \
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
