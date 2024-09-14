#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

if (NOT DEFINED Kokkos_CXX_COMPILER)
  message(FATAL_ERROR "ADIOS: Kokkos module requires the Kokkos_CXX_COMPILER variable")
endif()

# CXX Compiler settings only in for this subdir
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_COMPILER "${Kokkos_CXX_COMPILER}")

add_library(adios2_core_kokkos adiosKokkos.h adiosKokkos.cpp)

set_target_properties(adios2_core_kokkos PROPERTIES
  VISIBILITY_INLINES_HIDDEN ON
  INCLUDE_DIRECTORIES "$<BUILD_INTERFACE:${ADIOS2_SOURCE_DIR}/source>;$<BUILD_INTERFACE:${ADIOS2_BINARY_DIR}/source>"
  EXPORT_NAME core_kokkos
  OUTPUT_NAME adios2${ADIOS2_LIBRARY_SUFFIX}_core_kokkos
  )

kokkos_compilation(SOURCE adiosKokkos.cpp)
target_link_libraries(adios2_core_kokkos PRIVATE Kokkos::kokkos)
