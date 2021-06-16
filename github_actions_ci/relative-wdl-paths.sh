#!/bin/bash
set -e -o pipefail

echo Replacing relative paths to same directory
sed -i -- 's|import \"../tasks/|import \"|g' pipes/WDL/workflows/*.wdl
