#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="pugixml"
readonly ownership="pugixml Upstream <robot@adios2>"
readonly subtree="thirdparty/pugixml/pugixml"
readonly repo="https://github.com/zeux/pugixml.git"
readonly tag="master"
readonly shortlog="true"
readonly paths="
src/pugiconfig.hpp
src/pugixml.hpp
src/pugixml.cpp
"

extract_source () {
    git_archive
}

. "${BASH_SOURCE%/*}/../update-common.sh"
