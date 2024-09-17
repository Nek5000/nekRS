#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

# These tests should be *very* fast
set(CTEST_TEST_TIMEOUT 10)

set(BP3_DIR ${CMAKE_CURRENT_BINARY_DIR}/bp3)
set(BP4_DIR ${CMAKE_CURRENT_BINARY_DIR}/bp4)
set(BP5_DIR ${CMAKE_CURRENT_BINARY_DIR}/bp5)
set(BPfile_DIR ${CMAKE_CURRENT_BINARY_DIR}/bpfile)
set(FS_DIR ${CMAKE_CURRENT_BINARY_DIR}/filestream)
file(MAKE_DIRECTORY ${BP3_DIR})
file(MAKE_DIRECTORY ${BP4_DIR})
file(MAKE_DIRECTORY ${BP5_DIR})
file(MAKE_DIRECTORY ${BPfile_DIR})
file(MAKE_DIRECTORY ${FS_DIR})

gtest_add_tests_helper(Interface MPI_ALLOW ADIOS Interface.  .BP3
    WORKING_DIRECTORY ${BP3_DIR} EXTRA_ARGS "BP3")
gtest_add_tests_helper(Interface MPI_ALLOW ADIOS Interface.  .BPfile
    WORKING_DIRECTORY ${BPfile_DIR} EXTRA_ARGS "BPfile")
gtest_add_tests_helper(Write MPI_ALLOW ADIOSInterface Interface. "")
gtest_add_tests_helper(DefineVariable MPI_ALLOW ADIOS Interface. "")
gtest_add_tests_helper(DefineAttribute MPI_ALLOW ADIOS Interface. "")
gtest_add_tests_helper(Selection MPI_NONE ADIOS Interface.  .BP3
    WORKING_DIRECTORY ${BP3_DIR} EXTRA_ARGS "BP3")
gtest_add_tests_helper(Selection MPI_NONE ADIOS Interface.  .BPfile
    WORKING_DIRECTORY ${BPfile_DIR} EXTRA_ARGS "BPfile")
gtest_add_tests_helper(NoMpi MPI_NONE ADIOS Interface. "")
gtest_add_tests_helper(Info MPI_NONE ADIOS Interface. "")
