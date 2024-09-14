#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.6)
project(adios_plugin_test CXX)
enable_testing()

find_package(adios2 REQUIRED)

option(BUILD_SHARED_LIBS "build shared libs" ON)

#---------- Engine Plugin Tests

include(GenerateExportHeader)

add_library(PluginEngineWrite
  ../../../examples/plugins/engine/ExampleWritePlugin.cpp
)
target_link_libraries(PluginEngineWrite adios2::cxx11 adios2::core)
generate_export_header(PluginEngineWrite BASE_NAME plugin_engine_write)
target_include_directories(PluginEngineWrite PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
)

add_library(PluginEngineRead
  ../../../examples/plugins/engine/ExampleReadPlugin.cpp
)
target_link_libraries(PluginEngineRead adios2::cxx11 adios2::core)
generate_export_header(PluginEngineRead BASE_NAME plugin_engine_read)
target_include_directories(PluginEngineRead PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
)

if (WIN32)
  file(TO_NATIVE_PATH "${CMAKE_CURRENT_BINARY_DIR}/Debug" plugin_path)
else()
  set(plugin_path "${CMAKE_CURRENT_BINARY_DIR}")
endif()

# add write test
add_executable(adios_plugin_engine_write_test
  ../../../examples/plugins/engine/examplePluginEngineWrite.cpp
)
target_link_libraries(adios_plugin_engine_write_test adios2::cxx11)
add_test(NAME adios_plugin_engine_write_test COMMAND adios_plugin_engine_write_test)
set_tests_properties(adios_plugin_engine_write_test PROPERTIES
  ENVIRONMENT ADIOS2_PLUGIN_PATH=${plugin_path}
  )

# add read test
add_executable(adios_plugin_engine_read_test
  ../../../examples/plugins/engine/examplePluginEngineRead.cpp
)
target_link_libraries(adios_plugin_engine_read_test adios2::cxx11)
add_test(NAME adios_plugin_engine_read_test COMMAND adios_plugin_engine_read_test)
set_tests_properties(adios_plugin_engine_read_test PROPERTIES
  ENVIRONMENT ADIOS2_PLUGIN_PATH=${plugin_path}
  DEPENDS adios_plugin_engine_write_test
)

