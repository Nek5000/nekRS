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

cmake \
  -DCMAKE_INSTALL_PREFIX="${install_dir}" \
  -DBUILD_MPI=yes \
  -DBUILD_EXE=yes \
  -DBUILD_LIB=no \
  -DBUILD_DOC=no \
  -DLAMMPS_SIZES=smallbig \
  -DPKG_ADIOS=yes \
  "${source_dir}/cmake"
