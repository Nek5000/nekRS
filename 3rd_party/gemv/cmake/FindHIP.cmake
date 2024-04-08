# The MIT License (MIT)

# Copyright (c) 2014-2022 David Medina and Tim Warburton

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

if (NOT UNIX)
  message(FATAL_ERROR "HIP is only supported on Linux !")
endif()

# Search for HIP installation
if (NOT HIP_ROOT_DIR)
  # Search in user specified path first
  find_path(HIP_ROOT_DIR NAMES hipconfig PATHS ENV ROCM_PATH ENV HIP_PATH
    PATH_SUFFIXES bin DOC "HIP installed location" NO_DEFAULT_PATH)

  # Now search in default path
  find_path(HIP_ROOT_DIR NAMES hipconfig PATHS /opt/rocm PATH_SUFFIXES bin
    DOC "HIP installed location")

  # Check if we found HIP installation
  if(HIP_ROOT_DIR)
    # If so, fix the path
    string(REGEX REPLACE "[/\\\\]?bin[64]*[/\\\\]?$" "" HIP_ROOT_DIR
      ${HIP_ROOT_DIR})
    # And push it back to the cache
    set(HIP_ROOT_DIR ${HIP_ROOT_DIR} CACHE PATH "HIP installed location" FORCE)
  endif()
endif()

find_program(HIPCONFIG_EXECUTABLE NAMES hipconfig PATHS "${HIP_ROOT_DIR}"
  ENV ROCM_PATH ENV HIP_PATH /opt/rocm
  PATH_SUFFIXES bin NO_DEFAULT_PATH)

if (NOT HIPCONFIG_EXECUTABLE)
  find_program(HIPCONFIG_EXECUTABLE hipconfig)
endif()
mark_as_advanced(HIPCONFIG_EXECUTABLE)

if (HIPCONFIG_EXECUTABLE AND NOT HIP_PLATFORM)
  execute_process(COMMAND ${HIPCONFIG_EXECUTABLE} --platform
    OUTPUT_VARIABLE _hip_platform
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(HIP_PLATFORM ${_hip_platform} CACHE STRING
    "HIP platform as computed by hipconfig")
  mark_as_advanced(HIP_PLATFORM)
endif()

if (HIP_PLATFORM)
  if (${HIP_PLATFORM} STREQUAL "hcc" OR ${HIP_PLATFORM} STREQUAL "amd")
    set(HIP_LIBRARIES
      "${HIP_ROOT_DIR}/lib/libamdhip64.so;${HIP_ROOT_DIR}/lib/libhiprtc.so")
    set(HIP_RUNTIME_DEFINE "__HIP_PLATFORM_AMD__")
    set(HIP_INCLUDE_DIRS "${HIP_ROOT_DIR}/include")
  elseif (${HIP_PLATFORM} STREQUAL "nvidia")
    find_package(CUDA REQUIRED)
    find_package(CUDAToolkit REQUIRED)
    set(HIP_LIBRARIES "CUDA::cudart;CUDA::nvrtc")
      set(HIP_RUNTIME_DEFINE "__HIP_PLATFORM_NVIDIA__")
      set(HIP_INCLUDE_DIRS "${HIP_ROOT_DIR}/include;${CUDA_INCLUDE_DIRS}")
  endif()
endif()
mark_as_advanced(HIP_LIBRARIES)
mark_as_advanced(HIP_RUNTIME_DEFINE)
mark_as_advanced(HIP_INCLUDE_DIRS)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    HIP
    REQUIRED_VARS
    HIP_ROOT_DIR
    HIP_LIBRARIES
    HIP_RUNTIME_DEFINE
    HIP_INCLUDE_DIRS
    HIPCONFIG_EXECUTABLE
    HIP_PLATFORM)

if (NOT DEFINED ROCM_PATH)
  if (NOT DEFINED ENV{ROCM_PATH})
    set(ROCM_PATH "/opt/rocm" CACHE PATH "ROCm path")
  else()
    set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "ROCm path")
  endif()
endif()
set(CMAKE_MODULE_PATH "${ROCM_PATH}/lib/cmake" ${CMAKE_MODULE_PATH})

if (HIP_FOUND AND NOT TARGET gemv::HIP)
  add_library(gemv::HIP INTERFACE IMPORTED)
  set_target_properties(gemv::HIP PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "${HIP_RUNTIME_DEFINE}"
    INTERFACE_INCLUDE_DIRECTORIES "${HIP_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${HIP_LIBRARIES}"
  )

  find_package(hipblas)
  if (TARGET roc::hipblas)
    target_link_libraries(gemv::HIP INTERFACE roc::hipblas)
  endif()
endif()

if (TARGET roc::hipblas AND "${GEMV_DEFAULT_BACKEND}" STREQUAL "")
  set(GEMV_DEFAULT_BACKEND "hipblas")
endif()

if (TARGET gemv::HIP AND "${GEMV_DEFAULT_BACKEND}" STREQUAL "")
  set(GEMV_DEFAULT_BACKEND "hip")
endif()
