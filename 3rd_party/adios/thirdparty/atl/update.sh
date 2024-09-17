#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="atl"
readonly ownership="atl Upstream <robot@adios2>"
readonly subtree="thirdparty/atl/atl"
readonly repo="https://github.com/GTkorvo/atl.git"
readonly tag="master"
readonly shortlog="true"
readonly paths="
"

extract_source () {
    git_archive
}

. "${BASH_SOURCE%/*}/../update-common.sh"
