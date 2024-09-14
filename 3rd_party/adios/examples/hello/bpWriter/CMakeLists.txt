#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2HelloBPWriterExample)

if(NOT TARGET adios2_core)
  set(_components C)

  include(CheckLanguage)
  check_language(Fortran)
  if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
  endif()
  if(CMAKE_Fortran_COMPILER_LOADED)
    list(APPEND _components Fortran)
  endif()

  find_package(MPI COMPONENTS ${_components})
  if(MPI_FOUND)
    # Workaround for various MPI implementations forcing the link of C++ bindings
    add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

    list(APPEND _components MPI)
  endif()
  list(APPEND _components CXX)

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

add_executable(adios2_hello_bpWriter bpWriter.cpp)
target_link_libraries(adios2_hello_bpWriter adios2::cxx11)
install(TARGETS adios2_hello_bpWriter RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(adios2_hello_bpWriter_c bpWriter.c)
target_link_libraries(adios2_hello_bpWriter_c adios2::c)
install(TARGETS adios2_hello_bpWriter_c RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(adios2_hello_bpPutDeferred bpPutDeferred.cpp)
target_link_libraries(adios2_hello_bpPutDeferred adios2::cxx11)
install(TARGETS adios2_hello_bpPutDeferred RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(adios2_hello_bpSubStreams bpSubStreams.cpp)
target_link_libraries(adios2_hello_bpSubStreams adios2::cxx11)
install(TARGETS adios2_hello_bpSubStreams RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(ADIOS2_HAVE_Fortran)
  add_executable(adios2_hello_bpWriter_f bpWriter.F90)
  target_link_libraries(adios2_hello_bpWriter_f adios2::fortran)
  install(TARGETS adios2_hello_bpWriter_f RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

if(ADIOS2_HAVE_MPI)
  add_executable(adios2_hello_bpWriter_mpi bpWriter.cpp)
  target_link_libraries(adios2_hello_bpWriter_mpi adios2::cxx11_mpi MPI::MPI_C)
  install(TARGETS adios2_hello_bpWriter_mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(adios2_hello_bpWriter_c_mpi bpWriter.c)
  target_link_libraries(adios2_hello_bpWriter_c_mpi adios2::c_mpi MPI::MPI_C)
  install(TARGETS adios2_hello_bpWriter_c_mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(adios2_hello_bpPutDeferred_mpi bpPutDeferred.cpp)
  target_link_libraries(adios2_hello_bpPutDeferred_mpi adios2::cxx11_mpi MPI::MPI_C)
  install(TARGETS adios2_hello_bpPutDeferred_mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(adios2_hello_bpSubStreams_mpi bpSubStreams.cpp)
  target_link_libraries(adios2_hello_bpSubStreams_mpi adios2::cxx11_mpi MPI::MPI_C)
  install(TARGETS adios2_hello_bpSubStreams_mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  if(ADIOS2_HAVE_Fortran)
    add_executable(adios2_hello_bpWriter_f_mpi bpWriter.F90)
    target_link_libraries(adios2_hello_bpWriter_f_mpi adios2::fortran_mpi MPI::MPI_Fortran)
    install(TARGETS adios2_hello_bpWriter_f_mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
endif()
