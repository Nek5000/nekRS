#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.6)
project(adios_operator_plugin_test CXX)
enable_testing()

find_package(adios2 REQUIRED)

option(BUILD_SHARED_LIBS "build shared libs" ON)

set(ENV{ADIOS2_PLUGIN_PATH} "${adios2_DIR}/../../")

#---------- Operator Plugin Tests

# add write test
add_executable(adios_plugin_operator_write_test
  ../../../examples/plugins/operator/examplePluginOperatorWrite.cpp
)
target_link_libraries(adios_plugin_operator_write_test adios2::cxx11)
add_test(NAME adios_plugin_operator_write_test COMMAND adios_plugin_operator_write_test)

# add read test
add_executable(adios_plugin_operator_read_test
  ../../../examples/plugins/operator/examplePluginOperatorRead.cpp
)
target_link_libraries(adios_plugin_operator_read_test adios2::cxx11)
add_test(NAME adios_plugin_operator_read_test COMMAND adios_plugin_operator_read_test)
set_tests_properties(adios_plugin_operator_read_test PROPERTIES
    DEPENDS adios_plugin_operator_write_test)
