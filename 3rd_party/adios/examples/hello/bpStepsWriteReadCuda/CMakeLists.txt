#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2HelloBPStepsWriteReadCudaExample)

if(NOT TARGET adios2_core)
  set(_components CXX)

  find_package(CUDAToolkit REQUIRED)
  if(CUDAToolkit_FOUND)
    enable_language(CUDA)
  endif()

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

if(ADIOS2_HAVE_CUDA OR ADIOS2_HAVE_Kokkos_CUDA)
  add_executable(adios2_hello_bpStepsWriteReadCuda bpStepsWriteReadCuda.cu)
  target_link_libraries(adios2_hello_bpStepsWriteReadCuda PUBLIC adios2::cxx11 CUDA::cudart)
  set_target_properties(adios2_hello_bpStepsWriteReadCuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  install(TARGETS adios2_hello_bpStepsWriteReadCuda RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
