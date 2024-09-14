#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

file(REMOVE_RECURSE "${ADIOS2_BINARY_DIR}/testing/install/make/${TEST_CASE}")
file(COPY "${ADIOS2_SOURCE_DIR}/testing/install/${TEST_CASE}"
  DESTINATION "${ADIOS2_BINARY_DIR}/testing/install/make"
  )
if(ADIOS2_HAVE_MPI)
  set(have_mpi 1)
else()
  set(have_mpi 0)
endif()
set(ENV{PATH} "${ADIOS2_BINARY_DIR}/testing/install/install/${CMAKE_INSTALL_BINDIR}:$ENV{PATH}")
if(CMAKE_OSX_SYSROOT)
  set(isysroot "-isysroot ${CMAKE_OSX_SYSROOT}")
endif()
string(REPLACE ";" " " MPIEXEC_EXTRA_FLAGS "${MPIEXEC_EXTRA_FLAGS}")
execute_process(
  COMMAND "${MAKE_COMMAND}"
    "CC=${CMAKE_C_COMPILER}"
    "CXX=${CMAKE_CXX_COMPILER}"
    "FC=${CMAKE_Fortran_COMPILER}"
    "CFLAGS=${CMAKE_C_FLAGS}"
    "CXXFLAGS=${CMAKE_CXX_FLAGS}"
    "FFLAGS=${CMAKE_Fortran_FLAGS}"
    "LDFLAGS=${CMAKE_EXE_LINKER_FLAGS}"
    "MPIEXEC=${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_FLAGS}"
    "ADIOS2_HAVE_MPI=${have_mpi}"
    "ADIOS2_EXECUTABLE_SUFFIX=${ADIOS2_EXECUTABLE_SUFFIX}"
    "ISYSROOT=${isysroot}"
  WORKING_DIRECTORY "${ADIOS2_BINARY_DIR}/testing/install/make/${TEST_CASE}"
  RESULT_VARIABLE result
  )
if(result)
  message(FATAL_ERROR "Result of test was ${result}, should be 0")
endif()
