#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2HelloBPStepsWriteReadKokkosExample)

# CXX Compiler settings only in for this example
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT TARGET adios2_core)
  set(_components CXX)

  find_package(Kokkos 3.7 QUIET)
  if(Kokkos_FOUND AND DEFINED Kokkos_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER "${Kokkos_CXX_COMPILER}")
  endif()

  include(CheckLanguage)
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
  endif()
  if(CMAKE_Fortran_COMPILER_LOADED)
    list(APPEND _components Fortran)
  endif()

  find_package(MPI QUIET COMPONENTS ${_components})
  if(MPI_FOUND)
    # Workaround for various MPI implementations forcing the link of C++ bindings
    add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

    list(APPEND _components MPI)
  endif()

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
else()
  if(DEFINED Kokkos_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER "${Kokkos_CXX_COMPILER}")
  endif()
endif()

# flcl package needed for the Kokkos Fortran example
find_package(flcl QUIET)

# C++ Kokkos example using default layouts and different memory spaces
add_executable(adios2_hello_bpStepsWriteReadKokkos bpStepsWriteReadKokkos.cpp)
kokkos_compilation(SOURCE bpStepsWriteReadKokkos.cpp)
if(ADIOS2_HAVE_MPI)
  target_link_libraries(adios2_hello_bpStepsWriteReadKokkos adios2::cxx11_mpi MPI::MPI_C Kokkos::kokkos)
else()
  target_link_libraries(adios2_hello_bpStepsWriteReadKokkos adios2::cxx11 Kokkos::kokkos)
endif()
install(TARGETS adios2_hello_bpStepsWriteReadKokkos RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# C++ Kokkos example using default memory space and layouts
add_executable(adios2_hello_bpWriteReadKokkosView bpWriteReadKokkosView.cpp)
kokkos_compilation(SOURCE bpStepsWriteReadKokkos.cpp)
if(ADIOS2_HAVE_MPI)
  target_link_libraries(adios2_hello_bpWriteReadKokkosView adios2::cxx11_mpi MPI::MPI_C Kokkos::kokkos)
else()
  target_link_libraries(adios2_hello_bpWriteReadKokkosView adios2::cxx11 Kokkos::kokkos)
endif()
install(TARGETS adios2_hello_bpWriteReadKokkosView RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Fortran Kokkos example using default layouts and different memory spaces
if(ADIOS2_HAVE_MPI AND ADIOS2_HAVE_Fortran AND flcl_FOUND)
  add_executable(adios2_hello_bpStepsWriteReadKokkos_f bpStepsWriteReadKokkos.F90 view-f.f90 view-cxx.cc)
  target_link_libraries(adios2_hello_bpStepsWriteReadKokkos_f adios2::fortran_mpi MPI::MPI_Fortran flcl::flcl)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
    target_link_options(adios2_hello_bpStepsWriteReadKokkos_f PRIVATE LINKER:-lxlf90_r)
  endif()
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set_target_properties(adios2_hello_bpStepsWriteReadKokkos_f PROPERTIES LINKER_LANGUAGE Fortran)
  endif()
  install(TARGETS adios2_hello_bpStepsWriteReadKokkos_f RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
