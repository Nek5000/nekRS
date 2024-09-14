#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2SimulationsLorenzOdeExample)

if(NOT TARGET adios2_core)
  set(_components CXX)

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

add_executable(adios2_simulations_lorenzWriter lorenzWriter.cpp)
target_link_libraries(adios2_simulations_lorenzWriter adios2::cxx11)
install(TARGETS adios2_simulations_lorenzWriter DESTINATION ${CMAKE_INSTALL_BINDIR})

if(CMAKE_CXX_STANDARD GREATER_EQUAL 17)
  add_executable(adios2_simulations_lorenzReader lorenzReader.cpp)
  target_link_libraries(adios2_simulations_lorenzReader adios2::cxx11)
  install(TARGETS adios2_simulations_lorenzReader DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
