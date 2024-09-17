#!/bin/bash --login
# shellcheck disable=SC1091
set -e

source scripts/ci/gitlab-ci/setup-vars.sh

readonly CTEST_SCRIPT=scripts/ci/cmake/ci-${CI_JOB_NAME}.cmake
if [ ! -f "$CTEST_SCRIPT" ]
then
  echo "[E] Variable files does not exits: $CTEST_SCRIPT"
  exit 1
fi

readonly STEP=$1
if [ -z "$STEP" ]
then
  echo "[E] No argument given: $*"
  exit 2
fi

declare -a CTEST_STEP_ARGS=("-Ddashboard_full=OFF")
case ${STEP} in
  update)    CTEST_STEP_ARGS+=("${CI_UPDATE_ARGS}") ;;
  configure) CTEST_STEP_ARGS+=("-Ddashboard_do_submit=OFF") ;;
  build)     CTEST_STEP_ARGS+=("-Ddashboard_do_submit=OFF") ;;
  test)      CTEST_STEP_ARGS+=("-Ddashboard_do_submit=OFF") ;;
  submit)    CTEST_STEP_ARGS+=("-Ddashboard_do_submit_only=ON" "-Ddashboard_do_configure=ON" "-Ddashboard_do_build=ON" "-Ddashboard_do_test=ON") ;;
esac
CTEST_STEP_ARGS+=("-Ddashboard_do_${STEP}=ON")

echo "**********CTest Begin**********"
echo "ctest -VV -S ${CTEST_SCRIPT} ${CTEST_STEP_ARGS[*]}"
ctest -VV -S "${CTEST_SCRIPT}" "${CTEST_STEP_ARGS[@]}"
RET=$?
echo "**********CTest End************"

# EC: 0-127 this script errors, 128-INF ctest errors
if [ $RET -ne 0 ]
then
  (( RET += 127 ))
fi
exit $RET
