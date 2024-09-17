#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.6)
project(adios_c_test C)
enable_testing()

find_package(adios2 REQUIRED)

# Test using serial C bindings through adios2::c target.
add_executable(adios_c_test main_nompi.c)
target_link_libraries(adios_c_test adios2::c)
add_test(NAME adios_c_test COMMAND adios_c_test)

if(ADIOS2_HAVE_MPI)
  find_package(MPI REQUIRED)
  set(src main_mpi.c)
  set(mpilib MPI::MPI_C)
  set(mpiexec ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_FLAGS})

  # Test using mpi C bindings through adios2::c_mpi target.
  add_executable(adios_c_mpi_test main_mpi.c)
  target_link_libraries(adios_c_mpi_test adios2::c_mpi)
  add_test(NAME adios_c_mpi_test COMMAND ${mpiexec} $<TARGET_FILE:adios_c_mpi_test>)
else()
  set(src main_nompi.c)
  set(mpilib)
  set(mpiexec)
endif()

# Test using C bindings through legacy adios2::adios2 target.
add_executable(adios_adios2c_test ${src})
target_link_libraries(adios_adios2c_test adios2::adios2 ${mpilib})
add_test(NAME adios_adios2c_test COMMAND ${mpiexec} $<TARGET_FILE:adios_adios2c_test>)
