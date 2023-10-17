#!/bin/bash

# use sed to replace version strings of docker images based on versions defined in txt file
#
# skip this replacement for any version string line with the comment "#skip-global-version-pin"
#
# requires $MODULE_VERSIONS to be set to point to a text file with equal-sign-separated values
# export MODULE_VERSIONS="./requirements-modules.txt" && ./github_actions_ci/check-wdl-runtimes.sh

printf "Updating docker image tags in WDL files with those in ${MODULE_VERSIONS}\n\n"

while IFS='=' read module version; do
  OLD_TAG=$module
  NEW_TAG="$module:$version"
  NEW_TAG_BOLD="$module:$(tput bold)$version$(tput sgr0)"
  printf "Replacing: %-14s \n with tag: %-14s \n\n" "$OLD_TAG" "$NEW_TAG_BOLD"

  sed -i '' "/^\(.*\)[[:space:]]*#skip-global-version-pin[[:space:]]*$/!s|$OLD_TAG[^\"\']*|$NEW_TAG|g" pipes/WDL/tasks/*.wdl
done < $MODULE_VERSIONS

printf "Replacements skipped for lines marked with '#skip-global-version-pin' \n\n"