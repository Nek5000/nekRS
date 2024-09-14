#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2HelloHelloWorldExample)

if(NOT TARGET adios2_core)
  set(_components C CXX)

  find_package(MPI COMPONENTS C)
  if(MPI_FOUND)
    # Workaround for various MPI implementations forcing the link of C++ bindings
    add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

    list(APPEND _components MPI)
  endif()

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

add_executable(adios2_hello_helloWorld_c hello-world.c)
target_link_libraries(adios2_hello_helloWorld_c adios2::c)
install(TARGETS adios2_hello_helloWorld_c RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(adios2_hello_helloWorld hello-world.cpp)
target_link_libraries(adios2_hello_helloWorld adios2::cxx11)
install(TARGETS adios2_hello_helloWorld RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_executable(adios2_hello_helloWorld_hl hello-world-hl.cpp)
target_link_libraries(adios2_hello_helloWorld_hl adios2::cxx11)
install(TARGETS adios2_hello_helloWorld_hl RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(ADIOS2_HAVE_MPI)
  add_executable(adios2_hello_helloWorld_c_mpi hello-world.c)
  target_link_libraries(adios2_hello_helloWorld_c_mpi adios2::c_mpi MPI::MPI_C)
  install(TARGETS adios2_hello_helloWorld_c_mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(adios2_hello_helloWorld_mpi hello-world.cpp)
  target_link_libraries(adios2_hello_helloWorld_mpi adios2::cxx11_mpi MPI::MPI_C)
  install(TARGETS adios2_hello_helloWorld_mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(adios2_hello_helloWorld_hl_mpi hello-world-hl.cpp)
  target_link_libraries(adios2_hello_helloWorld_hl_mpi adios2::cxx11_mpi MPI::MPI_C)
  install(TARGETS adios2_hello_helloWorld_hl_mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
