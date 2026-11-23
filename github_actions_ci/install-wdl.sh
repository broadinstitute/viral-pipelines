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

# Cromwell & womtool: pull the latest daily "develop" build from Docker Hub and copy the
# jars out of the image. This is the same artifact Terra runs, so CI tests against what
# Terra actually uses rather than a stale pinned GitHub release (Cromwell no longer cuts
# GitHub releases -- they only publish daily builds to Docker Hub). The daily build
# publishes only Docker images -- no standalone jars are hosted anywhere -- so we
# exfiltrate the jars here rather than running cromwell/womtool inside docker (which would
# break docker-in-docker task execution).
extract_jar_from_docker () {
	_image=$1   # e.g. broadinstitute/cromwell:develop
	_src=$2     # path inside image, e.g. /app/cromwell.jar
	_dest=$3    # output filename in cwd, e.g. cromwell.jar
	echo "Pulling $_image to extract $_dest"
	docker pull --quiet "$_image"
	# --user keeps host file ownership as the runner (not root); cp -L dereferences the
	# versioned symlink (docker cp would copy a broken symlink since the jar filename
	# embeds a daily-changing git hash); readlink logs the exact build to CI output.
	docker run --rm \
		--user "$(id -u):$(id -g)" \
		-v "$PWD:/out" \
		--entrypoint /bin/bash \
		"$_image" \
		-c "echo -n 'extracted '; readlink -f '$_src'; cp -L '$_src' '/out/$_dest'"
}

extract_jar_from_docker broadinstitute/cromwell:develop /app/cromwell.jar cromwell.jar
extract_jar_from_docker broadinstitute/womtool:develop /app/womtool.jar womtool.jar

fetch_jar_from_github dnanexus dxCompiler dxCompiler 2.15.0

TGZ=dx-toolkit-v0.311.0-ubuntu-20.04-amd64.tar.gz
echo "Fetching $TGZ"
wget --quiet https://dnanexus-sdk.s3.amazonaws.com/$TGZ
tar -xzpf $TGZ