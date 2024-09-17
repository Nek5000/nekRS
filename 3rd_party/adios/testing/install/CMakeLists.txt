#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

function(add_install_cmake_test name)
  add_test(NAME Install.CMake.${name} COMMAND "${CMAKE_COMMAND}"
    "-DADIOS2_SOURCE_DIR=${CMAKE_SOURCE_DIR}"
    "-DADIOS2_BINARY_DIR=${CMAKE_BINARY_DIR}"
    "-DCMAKE_GENERATOR=${CMAKE_GENERATOR}"
    "-DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}"
    "-DCMAKE_GENERATOR_TOOLSET=${CMAKE_GENERATOR_TOOLSET}"
    "-DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}"
    "-DTEST_CASE=${name}"
    "-DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}"
    "-DCMAKE_INSTALL_CMAKEDIR=${CMAKE_INSTALL_CMAKEDIR}"
    "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
    "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
    "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
    "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
    "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
    "-DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}"
    "-DBUILD_TYPE=$<CONFIG>"
    "-DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}"
    "-DMPIEXEC_EXTRA_FLAGS=${MPIEXEC_EXTRA_FLAGS}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/run_cmake.cmake"
    )
  set_tests_properties(Install.CMake.${name} PROPERTIES
    DEPENDS Install.Setup
    FIXTURES_REQUIRED Install
    )
endfunction()

function(add_install_make_test name)
  if(ADIOS2_HAVE_MPI)
    set(cc "${MPI_C_COMPILER}")
    set(cxx "${MPI_CXX_COMPILER}")
    set(fc "${MPI_Fortran_COMPILER}")
    set(mpiexec "${MPIEXEC_EXECUTABLE}")
    set(mpiexec_extra_flags "${MPIEXEC_EXTRA_FLAGS}")
  else()
    set(cc "${CMAKE_C_COMPILER}")
    set(cxx "${CMAKE_CXX_COMPILER}")
    set(fc "${CMAKE_Fortran_COMPILER}")
    set(mpiexec)
    set(mpiexec_extra_flags)
  endif()
  add_test(NAME Install.Make.${name} COMMAND "${CMAKE_COMMAND}"
    "-DADIOS2_SOURCE_DIR=${CMAKE_SOURCE_DIR}"
    "-DADIOS2_BINARY_DIR=${CMAKE_BINARY_DIR}"
    "-DADIOS2_HAVE_MPI=${ADIOS2_HAVE_MPI}"
    "-DADIOS2_EXECUTABLE_SUFFIX=${ADIOS2_EXECUTABLE_SUFFIX}"
    "-DMAKE_COMMAND=${MAKE_COMMAND}"
    "-DCMAKE_C_COMPILER=${cc}"
    "-DCMAKE_CXX_COMPILER=${cxx}"
    "-DCMAKE_Fortran_COMPILER=${fc}"
    "-DCMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
    "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
    "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
    "-DCMAKE_EXE_LINKER_FLAGS=${CMAKE_EXE_LINKER_FLAGS}"
    "-DTEST_CASE=${name}"
    "-DCMAKE_INSTALL_BINDIR=${CMAKE_INSTALL_BINDIR}"
    "-DMPIEXEC_EXECUTABLE=${mpiexec}"
    "-DMPIEXEC_EXTRA_FLAGS=${mpiexec_extra_flags}"
    "-DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/run_make.cmake"
    )
  set_tests_properties(Install.Make.${name} PROPERTIES
    DEPENDS Install.Setup
    FIXTURES_REQUIRED Install
    )
endfunction()

add_test(NAME Install.Setup COMMAND "${CMAKE_COMMAND}"
  "-DADIOS2_SOURCE_DIR=${CMAKE_SOURCE_DIR}"
  "-DADIOS2_BINARY_DIR=${CMAKE_BINARY_DIR}"
  "-DBUILD_TYPE=$<CONFIG>"
  -P "${CMAKE_CURRENT_SOURCE_DIR}/run_install.cmake"
  )
set_tests_properties(Install.Setup PROPERTIES
  FIXTURES_SETUP Install
  )

# Remove -Wall -Werror in this scope
if(CMAKE_C_FLAGS)
  string(REPLACE "-Wall"    " " CMAKE_C_FLAGS   ${CMAKE_C_FLAGS})
  string(REPLACE "-Werror"  " " CMAKE_C_FLAGS   ${CMAKE_C_FLAGS})
endif()

if(CMAKE_CXX_FLAGS)
  string(REPLACE "-Wall"    " " CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  string(REPLACE "-Werror"  " " CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
endif()

add_install_cmake_test(C)
add_install_cmake_test(CXX11)
if(ADIOS2_HAVE_Fortran)
  add_install_cmake_test(Fortran)
endif()

if(NOT WIN32)
  find_program(MAKE_COMMAND make)
  if(MAKE_COMMAND)
    add_install_make_test(C)
    add_install_make_test(CXX11)
    if(ADIOS2_HAVE_Fortran)
      add_install_make_test(Fortran)
    endif()
  endif()
endif()

add_install_cmake_test(EnginePlugin)
if(ADIOS2_HAVE_Sodium)
  add_install_cmake_test(EncryptionOperator)
endif()
if (ADIOS2_HAVE_Catalyst)
  add_install_cmake_test(CatalystEnginePlugin)
endif()
