#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2HelloThreadWriteExample)

if(NOT TARGET adios2_core)
  set(_components CXX)

  find_package(Threads QUIET)

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

if(Threads_FOUND)
  add_executable(adios2_hello_bpThreadWrite bpThreadWrite.cpp)
  target_link_libraries(adios2_hello_bpThreadWrite adios2::cxx11 ${CMAKE_THREAD_LIBS_INIT})
  install(TARGETS adios2_hello_bpThreadWrite RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
