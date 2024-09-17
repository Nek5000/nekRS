#!/bin/bash

echo
echo
echo "****************************************"
echo "  Installing DILL"
echo "****************************************"

mkdir dill
cd dill
git clone https://github.com/GTKorvo/dill.git source
mkdir build
cd build
cmake \
  -DCMAKE_BUILD_TYPE=$1 \
  -DBUILD_TESTING=OFF \
  -DCMAKE_INSTALL_PREFIX=${PWD}/../install \
  ../source
cmake --build . -j4 --config $1
cmake --install . --config $1
