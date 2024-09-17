#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2HelloSstWriterExample)

if(NOT TARGET adios2_core)
  set(_components CXX)

  find_package(MPI COMPONENTS C)
  if(MPI_FOUND)
    # Workaround for various MPI implementations forcing the link of C++ bindings
    add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

    list(APPEND _components MPI)
  endif()

  if (NOT WIN32)
    find_package(LIBFABRIC 1.6 QUIET)
  endif ()

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

if(ADIOS2_HAVE_SST)
  add_executable(adios2_hello_sstWriter sstWriter.cpp)
  if(ADIOS2_HAVE_MPI)
    target_link_libraries(adios2_hello_sstWriter adios2::cxx11_mpi MPI::MPI_C)
  else()
    target_link_libraries(adios2_hello_sstWriter adios2::cxx11)
  endif()
  install(TARGETS adios2_hello_sstWriter RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
