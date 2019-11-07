#!/bin/bash

set -e

pushd docs
make html && echo "Docs built successfully!" || echo "Docs did NOT build successfully."
popd
