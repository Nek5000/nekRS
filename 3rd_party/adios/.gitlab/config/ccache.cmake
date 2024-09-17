##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set(version 4.6.1)
set(arch x86_64)

if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
  set(sha256sum da1e1781bc1c4b019216fa16391af3e1daaee7e7f49a8ec9b0cdc8a1d05c50e2)
  set(base_url https://github.com/ccache/ccache/releases/download)
  set(platform linux)
  set(extension tar.xz)
elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
  set(sha256sum 3e36ba8c80fbf7f2b95fe0227b9dd1ca6143d721aab052caf0d5729769138059)
  set(full_url https://gitlab.kitware.com/utils/ci-utilities/-/package_files/534/download)
  set(filename ccache)
  set(extension tar.gz)
elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
  set(sha256sum a6c6311973aa3d2aae22424895f2f968e5d661be003b25f1bd854a5c0cd57563)
  set(base_url https://github.com/ccache/ccache/releases/download)
  set(platform windows)
  set(extension zip)
else()
  message(FATAL_ERROR "Unrecognized platform ${CMAKE_HOST_SYSTEM_NAME}")
endif()

if(NOT DEFINED filename)
  set(filename "ccache-${version}-${platform}-${arch}")
endif()

set(tarball "${filename}.${extension}")

if(NOT DEFINED full_url)
  set(full_url "${base_url}/v${version}/${tarball}")
endif()

file(DOWNLOAD
  "${full_url}" $ENV{CCACHE_INSTALL_DIR}/${tarball}
  EXPECTED_HASH SHA256=${sha256sum}
  SHOW_PROGRESS
  )

execute_process(
  COMMAND ${CMAKE_COMMAND} -E tar xf ${tarball}
  WORKING_DIRECTORY $ENV{CCACHE_INSTALL_DIR}
  RESULT_VARIABLE extract_results
  )

if(extract_results)
  message(FATAL_ERROR "Extracting `${tarball}` failed: ${extract_results}.")
endif()

file(RENAME $ENV{CCACHE_INSTALL_DIR}/${filename} $ENV{CCACHE_INSTALL_DIR}/ccache)
