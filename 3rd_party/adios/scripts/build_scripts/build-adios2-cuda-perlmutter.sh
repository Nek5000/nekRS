#!/bin/bash
# shellcheck disable=SC2191

module load cudatoolkit/11.7
module load gcc/11.2.0
module load cmake/3.24.3
module refresh

######## User Configurations ########
ADIOS2_HOME=$(pwd)
BUILD_DIR=${ADIOS2_HOME}/build-cuda-perlmutter
INSTALL_DIR=${ADIOS2_HOME}/install-cuda-perlmutter

num_build_procs=4

######## ADIOS2 ########
mkdir -p "${BUILD_DIR}"
rm -f "${BUILD_DIR}/CMakeCache.txt"
rm -rf "${BUILD_DIR}/CMakeFiles"

ARGS_ADIOS=(
    -D CMAKE_INSTALL_PREFIX="${INSTALL_DIR}"
    -D CMAKE_CXX_COMPILER=g++
    -D CMAKE_C_COMPILER=gcc

    -D ADIOS2_USE_CUDA=ON

    -D CMAKE_POSITION_INDEPENDENT_CODE=TRUE
    -D BUILD_SHARED_LIBS=ON
    -D ADIOS2_USE_Fortran=OFF
)
cmake "${ARGS_ADIOS[@]}" -S "${ADIOS2_HOME}" -B "${BUILD_DIR}"
cmake --build "${BUILD_DIR}" -j${num_build_procs}
cmake --install "${BUILD_DIR}"
