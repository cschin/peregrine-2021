#!/bin/bash
set -ex
echo $1
aws s3 cp $1/run.sh run.sh
bash run.sh
