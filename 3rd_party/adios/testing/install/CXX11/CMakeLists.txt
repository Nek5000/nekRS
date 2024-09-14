#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.6)
project(adios_cxx11_test CXX)
enable_testing()

find_package(adios2 REQUIRED)

# Test using serial C++11 bindings through adios2::cxx11 target.
add_executable(adios_cxx11_test main_nompi.cxx)
target_link_libraries(adios_cxx11_test adios2::cxx11)
add_test(NAME adios_cxx11_test COMMAND adios_cxx11_test)

if(ADIOS2_HAVE_MPI)
  # Avoid using MPI::MPI_CXX
  enable_language(C)

  find_package(MPI REQUIRED)
  set(mpilib MPI::MPI_C)
  set(src main_mpi.cxx)
  set(mpiexec ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_FLAGS})
  add_compile_definitions(OMPI_SKIP_MPICXX MPICH_SKIP_MPICXX)

  # Test using mpi C++11 bindings through adios2::cxx11_mpi target.
  add_executable(adios_cxx11_mpi_test main_mpi.cxx)
  target_link_libraries(adios_cxx11_mpi_test adios2::cxx11_mpi)
  add_test(NAME adios_cxx11_mpi_test COMMAND ${mpiexec} $<TARGET_FILE:adios_cxx11_mpi_test>)
else()
  set(src main_nompi.cxx)
  set(mpilib)
  set(mpiexec)
endif()

# Test using C++11 bindings through legacy adios2::adios2 target.
add_executable(adios_adios2cxx11_test ${src})
target_link_libraries(adios_adios2cxx11_test adios2::adios2 ${mpilib})
add_test(NAME adios_adios2cxx11_test COMMAND ${mpiexec} $<TARGET_FILE:adios_adios2cxx11_test>)
