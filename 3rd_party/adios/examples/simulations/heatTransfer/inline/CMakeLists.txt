#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2SimulationsHeatTransferInlineExample)

if(NOT TARGET adios2_core)
  set(_components C CXX)

  find_package(MPI COMPONENTS ${_components})
  if(MPI_FOUND)
    # Workaround for various MPI implementations forcing the link of C++ bindings
    add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

    list(APPEND _components MPI)
  endif()

  find_package(Threads QUIET)

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
endif()

if(ADIOS2_HAVE_MPI AND Threads_FOUND)
  add_executable(adios2_simulations_heatTransferInline
    main.cpp
    InlineIO.cpp
    ../write/HeatTransfer.cpp
    ../write/Settings.cpp
  )
  target_link_libraries(adios2_simulations_heatTransferInline adios2::cxx11_mpi MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT})
  install(TARGETS adios2_simulations_heatTransferInline RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
