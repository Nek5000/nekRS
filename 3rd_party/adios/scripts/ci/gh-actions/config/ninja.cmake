cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

set(version 1.11.0)

if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
  set(sha256sum 9726e730d5b8599f82654dc80265e64a10a8a817552c34153361ed0c017f9f02)
  set(platform linux)
elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
  set(sha256sum 21915277db59756bfc61f6f281c1f5e3897760b63776fd3d360f77dd7364137f)
  set(platform mac)
elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
  set(sha256sum d0ee3da143211aa447e750085876c9b9d7bcdd637ab5b2c5b41349c617f22f3b)
  set(platform win)
else()
  message(FATAL_ERROR "Unrecognized platform ${CMAKE_HOST_SYSTEM_NAME}")
endif()

set(tarball "ninja-${platform}.zip")

file(DOWNLOAD
    "https://github.com/ninja-build/ninja/releases/download/v${version}/${tarball}" $ENV{CI_ROOT_DIR}/.local/bin/${tarball}
  EXPECTED_HASH SHA256=${sha256sum}
  SHOW_PROGRESS
  )

execute_process(
  COMMAND ${CMAKE_COMMAND} -E tar xf ${tarball}
  WORKING_DIRECTORY $ENV{CI_ROOT_DIR}/.local/bin
  RESULT_VARIABLE extract_results
  )

if(extract_results)
  message(FATAL_ERROR "Extracting `${tarball}` failed: ${extract_results}.")
endif()
