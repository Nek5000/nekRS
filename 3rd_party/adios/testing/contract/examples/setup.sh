#!/bin/bash

source_dir="/opt/adios2/source/examples"
build_dir=$(readlink -f "${PWD}")/build
install_dir=$(readlink -f "${PWD}")/install

export source_dir
export build_dir
export install_dir

echo "source_dir  = \"${source_dir}\""
echo "build_dir   = \"${build_dir}\""
echo "install_dir = \"${install_dir}\""
