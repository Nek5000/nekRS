#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="GoogleTest"
readonly ownership="Google Test Upstream <robot@adios2>"
readonly subtree="thirdparty/GTest/googletest"
readonly repo="https://github.com/google/googletest.git"
readonly tag="main"
readonly shortlog="true"
readonly paths="
"

extract_source () {
    git_archive
}

. "${BASH_SOURCE%/*}/../update-common.sh"
