#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2SimulationsGrayScottKokkosExample)

# CXX Compiler settings only in for this example
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT TARGET adios2_core)
  set(_components CXX)

  find_package(MPI COMPONENTS C)
  if(MPI_FOUND)
    # Workaround for various MPI implementations forcing the link of C++ bindings
    add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

    list(APPEND _components MPI)
  endif()

  find_package(Kokkos 3.7 QUIET)
  if(Kokkos_FOUND AND DEFINED Kokkos_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER "${Kokkos_CXX_COMPILER}")
  endif()

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
else()
  if(DEFINED Kokkos_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER "${Kokkos_CXX_COMPILER}")
  endif()
endif()

if(ADIOS2_HAVE_MPI AND ADIOS2_HAVE_Kokkos)
  add_executable(adios2_simulations_gray-scott-kokkos
    main.cpp
    gray-scott.cpp
    settings.cpp
    writer.cpp
    restart.cpp
  )
  kokkos_compilation(SOURCE gray-scott.cpp)
  target_link_libraries(adios2_simulations_gray-scott-kokkos adios2::cxx11_mpi MPI::MPI_C Kokkos::kokkos)
  install(TARGETS adios2_simulations_gray-scott-kokkos RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
