#!/bin/bash
API_BASE="https://api.github.com/repos/ornladios/adios2"
USER=${STATUS_ROBOT_NAME}
TOKEN=${STATUS_ROBOT_KEY}
COMMIT=${CIRCLE_SHA1}
CDASH_STATUS_CONTEXT="cdash"
SOURCE_DIR="$(readlink -f "${CIRCLE_WORKING_DIRECTORY}"/source)"

build_status_body() {
  cat <<EOF
{
  "state": "success",
  "target_url": "https://open.cdash.org/index.php?compare1=61&filtercount=1&field1=revision&project=ADIOS&showfilters=0&limit=100&value1=${COMMIT}&showfeed=0",
  "description": "Build and test results available on CDash",
  "context": "${CDASH_STATUS_CONTEXT}"
}
EOF
}

PYTHON_SCRIPT="${SOURCE_DIR}/scripts/ci/circle/findStatus.py"
curl -u "${USER}:${TOKEN}" "${API_BASE}/commits/${COMMIT}/statuses" | python3 "${PYTHON_SCRIPT}" --context ${CDASH_STATUS_CONTEXT}
exit_status=$?
if [ "$exit_status" -ne 0 ]
then
  echo "Need to post a status for context ${CDASH_STATUS_CONTEXT}"
  postBody="$(build_status_body)"
  postUrl="${API_BASE}/statuses/${COMMIT}"
  curl -u "${USER}:${TOKEN}" "${postUrl}" -H "Content-Type: application/json" -H "Accept: application/vnd.github.v3+json" -d "${postBody}"
fi
