#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="pybind11"
readonly ownership="PyBind11 Upstream <robot@adios2>"
readonly subtree="thirdparty/pybind11/pybind11"
readonly repo="https://github.com/pybind/pybind11.git"
readonly tag="master"
readonly shortlog="true"
readonly paths="
  include/pybind11
"

extract_source () {
    git_archive
}

. "${BASH_SOURCE%/*}/../update-common.sh"
