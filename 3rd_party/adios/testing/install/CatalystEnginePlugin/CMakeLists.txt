#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.6)
project(adios_catalyst_test CXX)
enable_testing()

find_package(adios2 REQUIRED)

option(BUILD_SHARED_LIBS "build shared libs" ON)

message("adios2 dir: ${adios2_DIR}")
message("adios2 root: ${adios2_ROOT}")
message("adios2 lib dir: ${adios2_LIBRARY_DIRS}")

# ----------- Test the ParaView Catalyst Plugin
# will need to figure out how to get the path of the plugin
add_executable(adios_plugin_engine_catalyst_test
    testCatalystEngine.cpp
)
target_link_libraries(adios_plugin_engine_catalyst_test adios2::cxx11)
set(ENV{ADIOS2_PLUGIN_PATH} ${adios2_DIR}/../../)
configure_file(test-catalyst-config.xml ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
add_test(NAME adios_plugin_engine_catalyst_test COMMAND
    adios_plugin_engine_catalyst_test test-catalyst-config.xml)
