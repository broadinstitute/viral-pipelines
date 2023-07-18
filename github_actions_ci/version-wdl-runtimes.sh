#!/bin/bash

# use sed to replace version strings of docker images based on versions defined in txt file

# requires $MODULE_VERSIONS to be set to point to a text file with equal-sign-separated values
# export MODULE_VERSIONS="./requirements-modules.txt" && ./github_actions_ci/version-wdl-runtimes.sh

while IFS='=' read module version; do
  OLD_TAG=$module
  if ! grep -q "sha256" <<< "$version"; then
    echo "$module is specified using image tag"
    NEW_TAG="$module:$version"
  else
    echo "$module is specified using image build hash"
    NEW_TAG="$module@$version"
  fi
  echo Replacing $OLD_TAG with $NEW_TAG in all task WDL files
  sed -i "s|$OLD_TAG[^\"\']*|$NEW_TAG|g" pipes/WDL/tasks/*.wdl
  
done < $MODULE_VERSIONS