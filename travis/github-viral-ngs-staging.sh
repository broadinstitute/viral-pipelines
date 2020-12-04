#!/bin/bash
set -e -o pipefail

if [ -z "$TRAVIS_PULL_REQUEST_BRANCH" ]; then

	travis/check-wdl-runtimes.sh
	travis/flatten-wdls.sh > /dev/null

	VERSION=$(travis/list-docker-tags.sh | cut -f 2 -d ":" | tail -1); echo "version - $VERSION"

	eval $(ssh-agent)
	openssl aes-256-cbc \
		-K $encrypted_fb18189f5cc1_key \
		-iv $encrypted_fb18189f5cc1_iv \
		-in travis/github-deploy-id_rsa.enc \
		-out travis/github-deploy-id_rsa \
		-d
	chmod 400 travis/github-deploy-id_rsa
	ssh-add travis/github-deploy-id_rsa
	ssh-add -l -E md5
	git clone git@github.com:broadinstitute/viral-ngs-staging.git

	cd viral-ngs-staging
	if [ -z "$TRAVIS_TAG" ]; then
		git checkout -B $TRAVIS_BRANCH
	fi

	rm -rf *; cp -a ../pipes ../travis/github-staging/* .
	../travis/dockstoreyml.sh pipes/WDL/flattened/*.wdl > .dockstore.yml
	git add -A -f
	git diff-index --quiet HEAD || git commit -q -m "CI push github.com/broadinstitute/viral-pipelines:$VERSION"

	if [ "$TRAVIS_BRANCH" = "master" -o -n "$TRAVIS_TAG" ]; then
		# for dockstore, don't bother tagging every branch commit that is non-master -- just git push the branch instead
		git tag $VERSION
		git push origin --tags
	fi

	if [ -z "$TRAVIS_TAG" ]; then
		# if TRAVIS_TAG is set, skip this, since a separate Travis build is already pushing master branch anyway
		git push -f -u origin $TRAVIS_BRANCH
	fi

fi
