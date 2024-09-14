#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- #
#Distributed under the OSI - approved Apache License, Version 2.0. See
#accompanying file Copyright.txt for details.
#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- #

cmake_minimum_required(VERSION 3.12)
project(ADIOS2HelloSSTKokkosExample)

#CXX Compiler settings only in for this example
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (NOT TARGET adios2_core)
  set(_components CXX)

  find_package(Kokkos 3.7 QUIET)
  if (Kokkos_FOUND AND DEFINED Kokkos_CXX_COMPILER)
      set(CMAKE_CXX_COMPILER "${Kokkos_CXX_COMPILER}")
  endif()

  find_package(ADIOS2 REQUIRED COMPONENTS ${_components})
else()
  if (DEFINED Kokkos_CXX_COMPILER)
      set(CMAKE_CXX_COMPILER "${Kokkos_CXX_COMPILER}")
  endif()
endif()

if (ADIOS2_HAVE_Kokkos)
  add_executable(adios2_hello_sstWriterKokkos sstWriterKokkos.cpp)
  add_executable(adios2_hello_sstReaderKokkos sstReaderKokkos.cpp)
  kokkos_compilation(SOURCE sstWriterKokkos.cpp)
  kokkos_compilation(SOURCE sstReaderKokkos.cpp)
  target_link_libraries(adios2_hello_sstWriterKokkos adios2::cxx11 Kokkos::kokkos)
  install(TARGETS adios2_hello_sstWriterKokkos RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  target_link_libraries(adios2_hello_sstReaderKokkos adios2::cxx11 Kokkos::kokkos)
  install(TARGETS adios2_hello_sstReaderKokkos RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
