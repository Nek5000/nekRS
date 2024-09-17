#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="perfstubs"
readonly ownership="PerfStubs Upstream <kwrobot@kitware.com>"
readonly subtree="thirdparty/perfstubs/perfstubs"
readonly repo="https://github.com/khuck/perfstubs.git"
readonly tag="master"

readonly shortlog="true"
readonly paths="
LICENSE
README.md
perfstubs_api/README.md
perfstubs_api/config.h.in
perfstubs_api/timer.h
perfstubs_api/timer.c
perfstubs_api/tool.h
"

extract_source () {
    git_archive
}

. "${BASH_SOURCE%/*}/../update-common.sh"
