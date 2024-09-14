#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="mingw-w64"
readonly ownership="mingw-w64 Upstream <robot@adios2>"
readonly subtree="thirdparty/mingw-w64/mingw-w64"
readonly repo="https://github.com/mirror/mingw-w64.git"
readonly tag="master"
readonly shortlog="true"
readonly paths="
mingw-w64-crt/misc/getopt.c
mingw-w64-headers/crt/getopt.h
"

extract_source () {
    git_archive
}

. "${BASH_SOURCE%/*}/../update-common.sh"
