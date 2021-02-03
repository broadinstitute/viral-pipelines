#!/bin/bash
# BEFORE RUNNING THIS, please read the README

# set project, bucket, and service account
PROJECT="pgs-automation"
BUCKET="viral-sequencing-automation"
SERVICE_ACCOUNT="launch-full-viral-illumina-wor@pgs-automation.iam.gserviceaccount.com"

FUNCTION="launch_workflow"

gcloud alpha functions deploy "${FUNCTION}" \
  --entry-point="${FUNCTION}" \
  --env-vars-file="config.yaml" \
  --max-instances="1" \
  --memory="1024MB" \
  --project="${PROJECT}" \
  --runtime="python37" \
  --service-account="${SERVICE_ACCOUNT}" \
  --timeout="300" \
  --trigger-event="google.storage.object.finalize" \
  --trigger-resource="${BUCKET}" \
  -q
