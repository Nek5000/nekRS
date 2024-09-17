#!/bin/bash
# shellcheck disable=SC2191

module load gcc/10.2
module load cuda/11.5
module load cmake/3.23
module refresh

######## User Configurations ########
ADIOS2_HOME=$(pwd)
BUILD_DIR=${ADIOS2_HOME}/build-cuda-summit
INSTALL_DIR=${ADIOS2_HOME}/install-cuda-summit

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
)
cmake "${ARGS_ADIOS[@]}" -S "${ADIOS2_HOME}" -B "${BUILD_DIR}"
cmake --build "${BUILD_DIR}" -j${num_build_procs}
cmake --install "${BUILD_DIR}"
