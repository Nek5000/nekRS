#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.6)
project(adios_fortran_test Fortran C)
enable_testing()

find_package(adios2 COMPONENTS Fortran REQUIRED)

# Test using serial Fortran bindings through adios2::fortran target.
add_executable(adios_fortran_test main_nompi_check.F90)
target_link_libraries(adios_fortran_test adios2::fortran)
add_test(NAME adios_fortran_test COMMAND adios_fortran_test)

if(ADIOS2_HAVE_MPI)
  find_package(MPI REQUIRED)
  set(src main_mpi_check.F90)
  set(mpilib MPI::MPI_Fortran)
  set(mpiexec ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_FLAGS})

  # Test using mpi Fortran bindings through adios2::fortran_mpi target.
  add_executable(adios_fortran_mpi_test main_mpi_check.F90)
  target_link_libraries(adios_fortran_mpi_test adios2::fortran_mpi)
  add_test(NAME adios_fortran_mpi_test COMMAND ${mpiexec} $<TARGET_FILE:adios_fortran_mpi_test>)
else()
  set(src main_nompi_check.F90)
  set(mpilib)
  set(mpiexec)
endif()

# Test using Fortran bindings through legacy adios2::adios2_f target.
add_executable(adios_adios2_f_test ${src})
target_link_libraries(adios_adios2_f_test adios2::adios2_f ${mpilib})
add_test(NAME adios_adios2_f_test COMMAND ${mpiexec} $<TARGET_FILE:adios_adios2_f_test>)
