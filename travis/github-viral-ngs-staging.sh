#!/bin/bash
set -e -o pipefail

openssl aes-256-cbc \
	-K $encrypted_fb18189f5cc1_key \
	-iv $encrypted_fb18189f5cc1_iv \
	-in travis/github-deploy-id_rsa.enc \
	-out travis/github-deploy-id_rsa \
	-d
chmod 400 travis/github-deploy-id_rsa
ssh-add travis/github-deploy-id_rsa
git clone git@github.com:broadinstitute/viral-ngs-staging.git
