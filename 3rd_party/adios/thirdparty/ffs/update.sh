#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="ffs"
readonly ownership="ffs Upstream <robot@adios2>"
readonly subtree="thirdparty/ffs/ffs"
readonly repo="https://github.com/GTkorvo/ffs.git"
readonly tag="master"
readonly shortlog="true"
readonly paths="
"

extract_source () {
    git_archive
}

. "${BASH_SOURCE%/*}/../update-common.sh"
