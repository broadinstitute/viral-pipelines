#!/bin/bash
set -e -o pipefail

fetch_jar_from_github () {
	_github_org=$1
	_repo_name=$2
	_tool_name=$3
	_jar_version=$4
	_jar_fname="$_tool_name-$_jar_version.jar"
	echo "Fetching $_jar_fname"
	wget --quiet https://github.com/$_github_org/$_repo_name/releases/download/$_jar_version/$_jar_fname
	ln -s $_jar_fname $_tool_name.jar
}

fetch_jar_from_github broadinstitute cromwell womtool 92
fetch_jar_from_github broadinstitute cromwell cromwell 92
fetch_jar_from_github dnanexus dxCompiler dxCompiler 2.15.0

TGZ=dx-toolkit-v0.311.0-ubuntu-20.04-amd64.tar.gz
echo "Fetching $TGZ"
wget --quiet https://dnanexus-sdk.s3.amazonaws.com/$TGZ
tar -xzpf $TGZ