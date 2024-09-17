#!/bin/bash --login
set -e

# Strip the job name prefix
export CI_JOB_NAME="${CI_JOB_NAME#*:}"
export CI_BUILD_NAME="${CI_COMMIT_BRANCH#github/}_${CI_JOB_NAME}"
export CI_ROOT_DIR="${CI_PROJECT_DIR}/.."
export CI_SITE_NAME="${GITLAB_SITE}"
export CI_SOURCE_DIR="${CI_PROJECT_DIR}"

if [ -z "$DOWNSTREAM_COMMIT_SHA" ]
then
  export CI_COMMIT_REF="${CI_COMMIT_SHA}"
else
  export CI_COMMIT_REF="${DOWNSTREAM_COMMIT_SHA}"
fi

if [ -z "$DOWNSTREAM_BRANCH_REF" ]
then
  export CI_BRANCH_REF="${CI_COMMIT_REF_NAME}"
else
  export CI_BRANCH_REF="${DOWNSTREAM_BRANCH_REF}"
fi

# In OLCF Crusher we must fix the build directory in the yml.
if [ -z "$CI_BIN_DIR" ]
then
  export CI_BIN_DIR="${CI_ROOT_DIR}/${CI_BUILD_NAME}"
fi

# In OLCF Gitlab our PRs branches tip commit is not the head commit of the PR,
# it is instead the so called merged_commit_sha as described in the GitHub Rest
# API for pull requests. We need to report to the CDASH the original commit
# thus, we set it here using the CTEST_UPDATE_VERSION_OVERRIDE CMake variable
if [[ ${CI_BRANCH_REF} =~ ^pr[0-9]+_.*$ ]]
then
  # Original commit it is always its 2nd parent
  ci_original_sha=$(git rev-parse "${CI_COMMIT_REF}^2")
  export CI_ORIGINAL_SHA="$ci_original_sha"
  export CI_UPDATE_ARGS="-DCTEST_UPDATE_VERSION_OVERRIDE=${CI_ORIGINAL_SHA}"
else
  export CI_ORIGINAL_SHA="${CI_COMMIT_REF}"
fi
