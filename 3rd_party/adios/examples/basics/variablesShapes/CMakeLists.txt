#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2BasicsVariablesShapesExample)

if(NOT TARGET adios2_core)
  set(_components CXX)

  find_package(MPI COMPONENTS C)
  if(MPI_FOUND)
    # Workaround for various MPI implementations forcing the link of C++ bindings
    add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

    list(APPEND _components MPI)
  endif()

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

add_executable(adios2_basics_variablesShapes variablesShapes.cpp)
target_link_libraries(adios2_basics_variablesShapes adios2::cxx11)
install(TARGETS adios2_basics_variablesShapes RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(adios2_basics_variablesShapes_hl variablesShapes_hl.cpp)
target_link_libraries(adios2_basics_variablesShapes_hl adios2::cxx11)
install(TARGETS adios2_basics_variablesShapes_hl RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(ADIOS2_HAVE_MPI)
  add_executable(adios2_basics_variablesShapes_mpi variablesShapes.cpp)
  target_link_libraries(adios2_basics_variablesShapes_mpi adios2::cxx11_mpi MPI::MPI_C)
  install(TARGETS adios2_basics_variablesShapes_mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(adios2_basics_variablesShapes_hl_mpi variablesShapes_hl.cpp)
  target_link_libraries(adios2_basics_variablesShapes_hl_mpi adios2::cxx11_mpi MPI::MPI_C)
  install(TARGETS adios2_basics_variablesShapes_hl_mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
