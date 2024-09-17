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

export CC=mpicc
export CXX=mpic++
export FC=mpifort

cmake \
  -DCMAKE_INSTALL_PREFIX="${install_dir}" \
  -DFPHSA_NAME_MISMATCHED=ON \
  -DPIO_ENABLE_TESTS=ON \
  -DPIO_ENABLE_EXAMPLES=ON \
  -DWITH_NETCDF=OFF \
  -DWITH_PNETCDF=ON \
  -DPnetCDF_PATH="$(spack location -i parallel-netcdf)" \
  -DWITH_ADIOS2=ON \
  "${source_dir}"
