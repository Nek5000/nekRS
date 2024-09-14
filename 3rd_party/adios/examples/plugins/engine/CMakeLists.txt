#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2PluginsEngineExample)

if(NOT TARGET adios2_core)
  set(_components CXX)

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

include(GenerateExportHeader)

add_library(adios2_plugins_exampleWritePlugin ExampleWritePlugin.cpp)
target_link_libraries(adios2_plugins_exampleWritePlugin adios2::cxx11 adios2_core)
generate_export_header(adios2_plugins_exampleWritePlugin BASE_NAME plugin_engine_write)
target_include_directories(adios2_plugins_exampleWritePlugin PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS adios2_plugins_exampleWritePlugin RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_library(adios2_plugins_exampleReadPlugin ExampleReadPlugin.cpp)
target_link_libraries(adios2_plugins_exampleReadPlugin adios2::cxx11 adios2_core)
generate_export_header(adios2_plugins_exampleReadPlugin BASE_NAME plugin_engine_read)
target_include_directories(adios2_plugins_exampleReadPlugin PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS adios2_plugins_exampleReadPlugin RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(adios2_plugins_examplePluginEngineWrite examplePluginEngineWrite.cpp)
target_link_libraries(adios2_plugins_examplePluginEngineWrite adios2::cxx11)
install(TARGETS adios2_plugins_examplePluginEngineWrite RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(adios2_plugins_examplePluginEngineRead examplePluginEngineRead.cpp)
target_link_libraries(adios2_plugins_examplePluginEngineRead adios2::cxx11)
install(TARGETS adios2_plugins_examplePluginEngineRead RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
