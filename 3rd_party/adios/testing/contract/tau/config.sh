#!/bin/bash

set -x
set -e

# shellcheck disable=SC1091
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/setup.sh"

# Fail if is not set
source_dir="${source_dir:?}"
build_dir="${build_dir:?}"
install_dir="${install_dir:?}"

mkdir -p "${build_dir}"
cd "${build_dir}"

cmake -DCMAKE_INSTALL_PREFIX="${install_dir}" "${source_dir}"
