#!/bin/bash

# use sed to replace version strings of docker images based on versions defined in txt file
#
# skip this replacement for any version string line with the comment "#skip-global-version-pin"
#
# requires $MODULE_VERSIONS to be set to point to a text file with equal-sign-separated values
# export MODULE_VERSIONS="./requirements-modules.txt" && ./github_actions_ci/version-wdl-runtimes.sh

# determine the directory containing this script and set CONTAINING_DIR variable

CONTAINING_DIR="$(realpath $(dirname $( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )))"

export MODULE_VERSIONS="${CONTAINING_DIR}/requirements-modules.txt"

printf "Updating docker image tags in WDL files with those in ${MODULE_VERSIONS}\n\n"

while IFS='=' read module version; do
  OLD_TAG=$module
  NEW_TAG="$module:$version"
  NEW_TAG_BOLD="$module:$(tput bold)$version$(tput sgr0)"
  printf "Replacing: %-14s \n with tag: %-14s \n\n" "$OLD_TAG" "$NEW_TAG_BOLD"

  sed -i '' "/^\(.*\)[[:space:]]*#skip-global-version-pin[[:space:]]*$/!s|$OLD_TAG[^\"\']*|$NEW_TAG|g" pipes/WDL/tasks/*.wdl
done < $MODULE_VERSIONS

printf "Replacements skipped for lines marked with '#skip-global-version-pin' \n\n"