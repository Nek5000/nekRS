#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2PluginsOperatorExample)

if(NOT TARGET adios2_core)
  set(_components CXX)

  find_package(Sodium QUIET)

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

if(ADIOS2_HAVE_Sodium)
  add_executable(adios2_plugins_exampleOperatorPluginWrite examplePluginOperatorWrite.cpp)
  target_link_libraries(adios2_plugins_exampleOperatorPluginWrite adios2::cxx11)
  install(TARGETS adios2_plugins_exampleOperatorPluginWrite RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(adios2_plugins_exampleOperatorPluginRead examplePluginOperatorRead.cpp)
  target_link_libraries(adios2_plugins_exampleOperatorPluginRead adios2::cxx11)
  install(TARGETS adios2_plugins_exampleOperatorPluginRead RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()