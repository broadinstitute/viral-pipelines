#!/bin/bash

set -e

pushd docs

make html

build_exit_code=$?
if [ $build_exit_code -eq 1 ]; then
    echo "Docs built successfully"
else
    echo "Docs did NOT build successfully"
    exit $build_exit_code
fi
popd
