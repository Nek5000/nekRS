#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ADIOS2SimulationsHeatTransferWriteExample)

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

if (ADIOS2_HAVE_MPI AND Threads_FOUND)
  add_executable(adios2_simulations_heatTransferWrite
    main.cpp
    HeatTransfer.cpp
    Settings.cpp
    IO_adios2.cpp
  )
  target_link_libraries(adios2_simulations_heatTransferWrite adios2::cxx11_mpi MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT})
  install(TARGETS adios2_simulations_heatTransferWrite RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  add_executable(adios2_simulations_heatTransferWrite_ascii
    main.cpp
    HeatTransfer.cpp
    Settings.cpp
    IO_ascii.cpp
  )
  target_link_libraries(adios2_simulations_heatTransferWrite_ascii MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT})
  install(TARGETS adios2_simulations_heatTransferWrite_ascii RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  if(ADIOS2_HAVE_HDF5)
    add_executable(adios2_simulations_heatTransferWrite_hdf5
      main.cpp
      HeatTransfer.cpp
      Settings.cpp
      IO_hdf5_a.cpp
    )
    target_include_directories(adios2_simulations_heatTransferWrite_hdf5 PRIVATE ${HDF5_C_INCLUDE_DIRS})
    target_link_libraries(adios2_simulations_heatTransferWrite_hdf5
      ${HDF5_C_LIBRARIES} MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT})
    install(TARGETS adios2_simulations_heatTransferWrite_hdf5 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

    if(HDF5_IS_PARALLEL)
      add_executable(adios2_simulations_heatTransferWrite_ph5
        main.cpp
        HeatTransfer.cpp
        Settings.cpp
        IO_ph5.cpp
      )
      target_include_directories(adios2_simulations_heatTransferWrite_ph5 PRIVATE ${HDF5_C_INCLUDE_DIRS})
      target_link_libraries(adios2_simulations_heatTransferWrite_ph5
        ${HDF5_C_LIBRARIES} MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT})
      install(TARGETS adios2_simulations_heatTransferWrite_ph5 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()

    if(NOT (HDF5_VERSION VERSION_LESS 1.11))
      add_executable(adios2_simulations_heatTransferWrite_h5mixer
        main.cpp
        HeatTransfer.cpp
        Settings.cpp
        IO_h5mixer.cpp
      )
      target_link_libraries(adios2_simulations_heatTransferWrite_h5mixer
        adios2::cxx11_mpi MPI::MPI_C ${CMAKE_THREAD_LIBS_INIT})
      install(TARGETS adios2_simulations_heatTransferWrite_h5mixer RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()
  endif()
endif()
