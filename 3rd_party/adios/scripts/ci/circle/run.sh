#!/bin/bash --login

# shellcheck source=/dev/null
. /etc/profile.d/modules.sh

# Parse the branch name used by the PR
API_BASE="https://api.github.com/repos/ornladios/adios2"
REALBRANCH="${CIRCLE_BRANCH}"
if [ -n "${CIRCLE_PR_NUMBER}" ]
then
  APIURL="${API_BASE}/pulls/${CIRCLE_PR_NUMBER}"
  RESULT="$(curl -s "${APIURL}" | python3 -c "import sys, json; print(json.load(sys.stdin)['head']['ref'])" 2> /dev/null)"
  exit_status=$?
  if [ "$exit_status" -eq 0 ]
  then
    REALBRANCH="$(echo "${RESULT}" | tr '/' '-')"
  fi
fi

export CI_SITE_NAME="Circle CI"
if [ -n "${CIRCLE_PR_NUMBER}" ]
then
  export CI_BUILD_NAME="pr${CIRCLE_PR_NUMBER}_${REALBRANCH}_${CIRCLE_JOB}"
else
  export CI_BUILD_NAME="${CIRCLE_BRANCH}_${CIRCLE_JOB}"
fi
export CI_COMMIT_REF=${CIRCLE_SHA1}
export CI_ROOT_DIR="${PWD}"
export CI_SOURCE_DIR="${CI_ROOT_DIR}/source"
export CI_BIN_DIR="${CI_ROOT_DIR}/${CIRCLE_JOB}"


STEP=$1
CTEST_SCRIPT=${CI_SOURCE_DIR}/scripts/ci/cmake/ci-${CIRCLE_JOB}.cmake

# Update and Test steps enable an extra step
CTEST_STEP_ARGS=""
case ${STEP} in
  update) CTEST_STEP_ARGS="${CTEST_STEP_ARGS} -Ddashboard_do_checkout=ON" ;;
  test) CTEST_STEP_ARGS="${CTEST_STEP_ARGS} -Ddashboard_do_end=ON" ;;
esac
CTEST_STEP_ARGS="${CTEST_STEP_ARGS} -Ddashboard_do_${STEP}=ON"

if [ -x /opt/cmake/bin/ctest ]
then
  CTEST=/opt/cmake/bin/ctest
elif [ -s /Applications/CMake.app/Contents/bin/ctest ]
then
  CTEST=/Applications/CMake.app/Contents/bin/ctest
else
  CTEST=ctest
fi

# OpenMPI specific setup
if [[ "${SYSTEM_JOBNAME}" =~ openmpi ]]
then
  # Workaround to quiet some warnings from OpenMPI
  export OMPI_MCA_btl_base_warn_component_unused=0
  export OMPI_MCA_btl_vader_single_copy_mechanism=none

  # https://github.com/open-mpi/ompi/issues/6518
  export OMPI_MCA_btl=self,tcp

  # Enable overscription in OpenMPI
  export OMPI_MCA_rmaps_base_oversubscribe=1
  export OMPI_MCA_hwloc_base_binding_policy=none
fi

echo "**********Env Begin**********"
env | sort
echo "**********Env End************"

echo "**********CTest Begin**********"
${CTEST} --version
echo ${CTEST} -VV -S "${CTEST_SCRIPT}" -Ddashboard_full=OFF "${CTEST_STEP_ARGS}"
# shellcheck disable=SC2086
${CTEST} -VV -S "${CTEST_SCRIPT}" -Ddashboard_full=OFF ${CTEST_STEP_ARGS}
RET=$?
echo "**********CTest End************"

exit ${RET}
