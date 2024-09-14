#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2BasicsGlobalArray1DExample)

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

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

if(ADIOS2_HAVE_MPI)
  add_library(adios2_basics_globalArray1D_mpivars_c OBJECT mpivars.c)
  target_link_libraries(adios2_basics_globalArray1D_mpivars_c MPI::MPI_C)
  install(TARGETS adios2_basics_globalArray1D_mpivars_c RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_library(adios2_basics_globalArray1D_decomp_c OBJECT decomp.c)
  target_link_libraries(adios2_basics_globalArray1D_decomp_c
    adios2_basics_globalArray1D_mpivars_c MPI::MPI_C)
  install(TARGETS adios2_basics_globalArray1D_decomp_c RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(adios2_basics_globalArray1DWrite_c globalArray1DWrite.c)
  target_link_libraries(adios2_basics_globalArray1DWrite_c
    adios2_basics_globalArray1D_mpivars_c
    adios2_basics_globalArray1D_decomp_c
    adios2::c_mpi MPI::MPI_C)
  install(TARGETS adios2_basics_globalArray1DWrite_c RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(adios2_basics_globalArray1DRead_c globalArray1DRead.c)
  target_link_libraries(adios2_basics_globalArray1DRead_c
    adios2_basics_globalArray1D_mpivars_c
    adios2_basics_globalArray1D_decomp_c
    adios2::c_mpi MPI::MPI_C)
  install(TARGETS adios2_basics_globalArray1DRead_c RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  if(ADIOS2_HAVE_Fortran)
    add_library(adios2_basics_globalArray1D_mpivars_f OBJECT mpivars.F90)
    target_link_libraries(adios2_basics_globalArray1D_mpivars_f MPI::MPI_Fortran)
    install(TARGETS adios2_basics_globalArray1D_mpivars_f RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

    add_library(adios2_basics_globalArray1D_decomp_f OBJECT decomp.F90)
    target_link_libraries(adios2_basics_globalArray1D_decomp_f
      adios2_basics_globalArray1D_mpivars_f MPI::MPI_Fortran)
    install(TARGETS adios2_basics_globalArray1D_decomp_f RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

    add_executable(adios2_basics_globalArray1DWrite_f globalArray1DWrite.F90)
    target_link_libraries(adios2_basics_globalArray1DWrite_f
      adios2_basics_globalArray1D_mpivars_f
      adios2_basics_globalArray1D_decomp_f
      adios2::fortran_mpi MPI::MPI_Fortran)
    install(TARGETS adios2_basics_globalArray1DWrite_f RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

    add_executable(adios2_basics_globalArray1DRead_f globalArray1DRead.F90)
    target_link_libraries(adios2_basics_globalArray1DRead_f
      adios2_basics_globalArray1D_mpivars_f
      adios2_basics_globalArray1D_decomp_f
      adios2::fortran_mpi MPI::MPI_Fortran)
    install(TARGETS adios2_basics_globalArray1DRead_f RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
endif()
