#!/bin/bash

if [ $# -eq 0 ]
then
  echo "$0 org/repo [SHA [TOKEN]]"
  exit 1
fi

GH_REPO=$1
if [ $# -gt 1 ]
then
  GH_SHA=$2
  if [ $# -gt 2 ]
  then
    GH_TOKEN=$3
    if [ $# -gt 3 ]
    then
      echo "$0 org/repo [SHA [TOKEN]]"
      exit 1
    fi
  fi
fi

if [ -z "${GH_SHA}" ]
then
  echo "Error: GH_SHA is undefined"
  exit 2
fi
if [ -z "${GH_TOKEN}" ]
then
  echo "Error: GH_TOKEN is undefined"
  exit 3
fi

cdash_url() {
  printf 'https://open.cdash.org/index.php?project=%s&subproject=%s&filtercount=1&showfilters=1&field1=revision&compare1=61&value1=%s' $1 $2 $3
}

build_status_body() {
cat <<EOF
{
  "state": "success",
  "target_url": "$(cdash_url GTKorvo DILL ${GH_SHA})",
  "description": "Build and test results available on CDash",
  "context": "CDash"
}
EOF
}

build_status_body | \
curl \
  --request POST \
  --url https://api.github.com/repos/${GH_REPO}/statuses/${GH_SHA} \
  --header "Accept: application/vnd.github.v3+json" \
  --header "authorization: Bearer ${GH_TOKEN}" \
  -d @-
