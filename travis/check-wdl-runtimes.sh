#!/bin/bash
# requires $MODULE_VERSIONS to be set to point to a text file with equal-sign-separated values
# export MODULE_VERSIONS="./requirements-modules.txt" && ./github_actions_ci/check-wdl-runtimes.sh

echo "Checking wdl container versions against ${MODULE_VERSIONS}"

# this is the newer script that simply validates existing version strings
should_error=false
for task_file in $(ls -1 pipes/WDL/tasks/*.wdl); do
    echo "Checking ${task_file}"
    while IFS='=' read module version; do
        OLD_TAG=$module
        NEW_TAG="$module:$version"
        
        offending_lines="$(grep -nE "^[^#]*$OLD_TAG" "${task_file}" | grep -v $NEW_TAG)"

        # if the expected tag is not seen, let us know the file and exit
        if [ $? -eq 0 ]; then
           offending_lines="$(echo "${offending_lines}" | sed 's/^/      /g')"
           echo "  \"$NEW_TAG\" needed in \"${task_file}\":"
           echo "$offending_lines"
           should_error=true #exit 1 eventually, but only after printing all of the problems
        fi
        offending_lines=""
        
    done < "${MODULE_VERSIONS}"
done
if [ "$should_error" = true ]; then
    exit 1
fi