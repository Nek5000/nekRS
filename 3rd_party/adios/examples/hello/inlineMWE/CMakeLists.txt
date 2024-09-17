#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2HelloInlineMWEExample)

if(NOT TARGET adios2_core)
  set(_components CXX)

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

add_executable(adios2_hello_inlineMWE inlineMWE.cpp)
target_link_libraries(adios2_hello_inlineMWE adios2::cxx11)
install(TARGETS adios2_hello_inlineMWE RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
