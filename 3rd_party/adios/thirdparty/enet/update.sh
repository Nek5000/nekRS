#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="enet"
readonly ownership="enet Upstream <robot@adios2>"
readonly subtree="thirdparty/enet/enet"
readonly repo="https://github.com/GTkorvo/enet.git"
readonly tag="master"
readonly shortlog="true"
readonly paths="
"

extract_source () {
    git_archive
}

. "${BASH_SOURCE%/*}/../update-common.sh"
