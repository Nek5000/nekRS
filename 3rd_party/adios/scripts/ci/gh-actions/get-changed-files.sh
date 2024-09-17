#!/bin/bash

case "${GITHUB_EVENT_NAME}"
in
  pull_request)
    BASE_SHA=$(jq -r .pull_request.base.sha "${GITHUB_EVENT_PATH}")
    HEAD_SHA=$(jq -r .pull_request.head.sha "${GITHUB_EVENT_PATH}")
    ;;
  push)
    BASE_SHA=$(jq -r .before "${GITHUB_EVENT_PATH}")
    HEAD_SHA=$(jq -r .after "${GITHUB_EVENT_PATH}")
    ;;
  *)
    echo "Unable to get changed files from '${GITHUB_EVENT_NAME}' event"
    exit 1
esac

echo "Event: ${GITHUB_EVENT_NAME}"
echo "Base: ${BASE_SHA}"
echo "Head: ${HEAD_SHA}"
echo ""

git fetch origin "${BASE_SHA}"

echo ""
echo "::group::All changed files"
git diff --name-only "${BASE_SHA}"..."${HEAD_SHA}" | tee all-changed-files.txt

echo "::group::Filtered changes"
grep -v '^docs/' all-changed-files.txt | tee filtered-changed-files.txt

echo "::group::Ignored changes"
grep '^docs/' all-changed-files.txt | tee ignored-changed-files.txt
