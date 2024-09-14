#!/bin/bash

readlinkf() { perl -MCwd -e 'print Cwd::abs_path shift' "$1"; }

set -e

STEP=$1
shift

export CI_SITE_NAME="GitHub Actions"
if [ "${GITHUB_EVENT_NAME}" = "pull_request" ]
then
  GH_PR_NUMBER=$(expr "${GITHUB_REF}" : 'refs/pull/\([^/]*\)')
  export CI_BUILD_NAME="pr${GH_PR_NUMBER}_${GITHUB_HEAD_REF}_${GH_YML_JOBNAME}_${GH_YML_BUILDTYPE}"
else
  export CI_BUILD_NAME="${GITHUB_REF#refs/heads/}_${GH_YML_JOBNAME}_${GH_YML_BUILDTYPE}"
fi
if [[ "${RUNNER_OS}" =~ "Windows" ]]
then
  export CI_ROOT_DIR="${GITHUB_WORKSPACE//\\//}"
  export CI_SOURCE_DIR="${GITHUB_WORKSPACE//\\//}/source"
else
  export CI_ROOT_DIR="${GITHUB_WORKSPACE}"
  export CI_SOURCE_DIR="${GITHUB_WORKSPACE}/source"
fi
export CI_BIN_DIR="${CI_ROOT_DIR}/build"

export CI_COMMIT_SHA=${GH_YML_SHA}

if command -v ctest3 >/dev/null
then
  CTEST=ctest3
else
  CTEST=ctest
fi

# Update and Test steps enable an extra step
CTEST_STEP_ARGS=""
case ${STEP} in
  test) CTEST_STEP_ARGS="${CTEST_STEP_ARGS} -Ddashboard_do_end=ON" ;;
esac
CTEST_STEP_ARGS="${CTEST_STEP_ARGS} -Ddashboard_do_${STEP}=ON"

echo "**********  Environment  **********"

env | sort

echo "********** Running CTest **********"

${CTEST} -VV \
  -S ${CI_SOURCE_DIR}/scripts/ci/cmake/${GH_YML_JOBNAME}.cmake \
  -DCTEST_BUILD_CONFIGURATION=${GH_YML_BUILDTYPE} \
  -Ddashboard_do_submit=ON \
  -Ddashboard_full=OFF \
  ${CTEST_STEP_ARGS} \
  "$@"

