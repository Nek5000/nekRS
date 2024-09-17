#!/usr/bin/env bash

set -e
set -x
shopt -s dotglob

readonly name="KWSys"
readonly ownership="KWSys Upstream <kwrobot@kitware.com>"
readonly subtree="thirdparty/KWSys/adios2sys"
readonly repo="https://gitlab.kitware.com/utils/kwsys.git"
readonly tag="master"

readonly shortlog="true"
readonly paths="
"

git_archive_no_attributes () {
    find . -name gitattributes | xargs git rm -f
    git archive --prefix="$name-reduced/" HEAD -- $paths | \
        tar -C "$extractdir" -x
}

extract_source () {
    git_archive_no_attributes
}

. "${BASH_SOURCE%/*}/../update-common.sh"
