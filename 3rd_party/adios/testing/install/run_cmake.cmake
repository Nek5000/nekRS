#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

file(REMOVE_RECURSE "${ADIOS2_BINARY_DIR}/testing/install/cmake/${TEST_CASE}")
if(WIN32)
  set(ENV{PATH} "${ADIOS2_BINARY_DIR}/testing/install/install/${CMAKE_INSTALL_BINDIR};$ENV{PATH}")
endif()
execute_process(COMMAND "${CMAKE_CTEST_COMMAND}"
  --build-and-test
  "${ADIOS2_SOURCE_DIR}/testing/install/${TEST_CASE}"
  "${ADIOS2_BINARY_DIR}/testing/install/cmake/${TEST_CASE}"
  --build-generator "${CMAKE_GENERATOR}"
  --build-generator-platform "${CMAKE_GENERATOR_PLATFORM}"
  --build-generator-toolset "${CMAKE_GENERATOR_TOOLSET}"
  --build-makeprogram "${CMAKE_MAKE_PROGRAM}"
  --build-options
    "-Dadios2_DIR=${ADIOS2_BINARY_DIR}/testing/install/install/${CMAKE_INSTALL_CMAKEDIR}"
    "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
    "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
    "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
    "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
    "-DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}"
    "-DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}"
    "-DMPIEXEC_EXTRA_FLAGS=${MPIEXEC_EXTRA_FLAGS}"
  --test-command "${CMAKE_CTEST_COMMAND}" -V
  -C "${BUILD_TYPE}"
  RESULT_VARIABLE result
  )
if(result)
  message(FATAL_ERROR "Result of test was ${result}, should be 0")
endif()
