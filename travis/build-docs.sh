#!/bin/bash

set -e

if [[ ( -n $TRAVIS_PULL_REQUEST && $TRAVIS_PULL_REQUEST != "false" ) || $TRAVIS_BRANCH = "master" || -n "$TRAVIS_TAG" ]]; then

	while IFS='=' read module version; do
		git clone --depth=1 https://github.com/$module.git
		pushd $module
		git checkout -qf $version
		popd
		cp -a $module/* .
	done < $MODULE_VERSIONS

    pushd docs
    make html && echo "Docs built successfully!" || echo "Docs did NOT build successfully."
    popd
fi