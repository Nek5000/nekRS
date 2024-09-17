#!/bin/bash

set -x
set -e

# shellcheck disable=SC1091
source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/setup.sh"

# Fail if is not set
build_dir="${build_dir:?}"

cmake --install "${build_dir}"
