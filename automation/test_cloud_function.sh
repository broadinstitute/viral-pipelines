#!/bin/bash

if [ $# -lt 1 ]; then
    echo "USAGE: ./test_cloud_function.sh BUCKET"
    exit 1
fi


BUCKET=$1

# add a file to the bucket
echo "hello world!" > test.txt
gsutil cp test.txt gs://${BUCKET}

echo "now visit https://console.cloud.google.com/logs/ to view the logs"

# clean up
rm test.txt
