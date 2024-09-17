#!/bin/bash

if [ "${GITHUB_EVENT_NAME}" = "pull_request" ]
then
  if [ -z "${BASE_REF}" ]
  then
    BASE_REF="$(jq -r .pull_request.base.ref "${GITHUB_EVENT_PATH}")"
  fi
  echo "Base ref: ${BASE_REF}"
  echo "Head ref: ${GITHUB_HEAD_REF}"

  if [ "${GITHUB_HEAD_REF}" = "${BASE_REF}" ]
  then
    echo "ERROR: Head ref must be distinctly different from the target branch"
    exit 1
  elif [[ ${GITHUB_HEAD_REF} =~ ^(master|release|release_[0-9]+)$ ]]
  then
    echo "ERROR: Head ref cannot be master, release, or release_XY"
    exit 2
  fi
  echo "OK"
else
  echo "Head ref name check skipped for ${GITHUB_EVENT_NAME} event"
fi
