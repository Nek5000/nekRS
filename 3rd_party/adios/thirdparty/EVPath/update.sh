#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="EVPath"
readonly ownership="EVPath Upstream <robot@adios2>"
readonly subtree="thirdparty/EVPath/EVPath"
readonly repo="https://github.com/GTkorvo/EVPath.git"
readonly tag="master"
readonly shortlog="true"
readonly paths="
"

extract_source () {
    git_archive
}

. "${BASH_SOURCE%/*}/../update-common.sh"
