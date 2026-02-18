#!/bin/bash
# requires $MODULE_VERSIONS to be set to point to a text file with equal-sign-separated values
# export MODULE_VERSIONS="./requirements-modules.txt" && ./github_actions_ci/check-wdl-runtimes.sh

function absolute_path() {
    local SOURCE="$1"
    while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
        DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
        if [[ "$OSTYPE" == "darwin"* ]]; then
            SOURCE="$(readlink "$SOURCE")"
        else
            SOURCE="$(readlink -f "$SOURCE")"
        fi
        [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
    done
    echo "$SOURCE"
}
SOURCE="${BASH_SOURCE[0]}"
SCRIPT=$(absolute_path "$SOURCE")
SCRIPT_DIRNAME="$(dirname "$SOURCE")"
SCRIPTPATH="$(cd -P "$(echo "$SCRIPT_DIRNAME")" &> /dev/null && pwd)"
SCRIPT="$SCRIPTPATH/$(basename "$SCRIPT")"
REPO_PATH="$(realpath "$SCRIPTPATH/../")"

# set MODULE_VERSIONS with default value of "$(realpath "$SCRIPTPATH/../requirements-modules.txt")", assumed to be located one level above this script
MODULE_VERSIONS="${MODULE_VERSIONS:-"${REPO_PATH}/requirements-modules.txt"}"

echo "Checking wdl container versions against ${MODULE_VERSIONS}"

# this is the newer script that simply validates existing version strings
should_error=false
for task_file in $(ls -1 pipes/WDL/tasks/*.wdl); do
    echo "Checking ${task_file}"
    while IFS='=' read module version; do
    	OLD_TAG=$module
    	NEW_TAG="$module:$version"

        # Escape dots in version for regex and build pattern that allows optional flavor suffix
        escaped_version=$(echo "$version" | sed 's/\./\\./g')
        # Pattern matches: module:version or module:version-flavor, followed by quote
        VALID_PATTERN="${module}:${escaped_version}(-[a-zA-Z0-9_]+)?[\"']"

        offending_lines="$(grep -nE "^[^#]*$OLD_TAG" "${task_file}" | grep -v '#skip-global-version-pin' | grep -Ev "$VALID_PATTERN")"

        # if the expected tag is not seen, let us know the file and exit
        if [ $? -eq 0 ]; then
           offending_lines="$(echo "${offending_lines}" | sed 's/^/      /g')"
           echo "  \"$NEW_TAG\" (or $NEW_TAG-<flavor>) needed in \"${task_file}\":"
           echo "$offending_lines"
           should_error=true #exit 1 eventually, but only after printing all of the problems
        fi
        offending_lines=""

    done < "${MODULE_VERSIONS}"
done
if [ "$should_error" = true ]; then
    exit 1
fi