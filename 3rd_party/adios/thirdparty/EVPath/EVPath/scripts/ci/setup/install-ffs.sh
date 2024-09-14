#!/bin/bash

echo
echo
echo "****************************************"
echo "  Installing FFS"
echo "****************************************"

export atl_ROOT="${PWD}/atl/install"

if [[ "${OS}" =~ "Windows" ]]
then
  extra_cmake_args="-DFFS_USE_DILL:BOOL=OFF"
else
  extra_cmake_args="-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON"
  export dill_ROOT="${PWD}/dill/install"
fi

mkdir ffs
cd ffs
git clone https://github.com/GTKorvo/ffs.git source
mkdir build
cd build
cmake ${extra_cmake_args} \
  -DCMAKE_BUILD_TYPE=$1 \
  -DBUILD_TESTING=OFF \
  -DCMAKE_INSTALL_PREFIX=${PWD}/../install \
  ../source
cmake --build . -j4 --config $1
cmake --install . --config $1
