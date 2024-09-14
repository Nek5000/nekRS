#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#
add_library(sst OBJECT
  dp/dp.c
  dp/dummy_dp.c
  dp/evpath_dp.c
  cp/cp_reader.c
  cp/cp_writer.c
  cp/cp_common.c
  cp/ffs_marshal.c
  sst_comm.cpp
)
if(BUILD_SHARED_LIBS)
  set_property(TARGET sst PROPERTY POSITION_INDEPENDENT_CODE 1)
endif()

target_link_libraries(sst PRIVATE adios2::thirdparty::perfstubs-interface)

if(ADIOS2_SST_HAVE_LIBFABRIC)
  target_sources(sst PRIVATE dp/rdma_dp.c)
  target_link_libraries(sst PRIVATE libfabric::libfabric)
  set(CMAKE_REQUIRED_INCLUDES ${LIBFABRIC_INCLUDE_DIRS})
  CHECK_INCLUDE_FILE(rdma/fi_ext_gni.h ADIOS2_SST_HAVE_FI_GNI)
  if(ADIOS2_SST_HAVE_CRAY_DRC)
    target_link_libraries(sst PRIVATE craydrc::craydrc)
  endif()
endif()

if(ADIOS2_SST_HAVE_UCX)
  target_sources(sst PRIVATE dp/ucx_dp.c)
  target_link_libraries(sst
                  PRIVATE ucx::ucx)
endif()

if(ADIOS2_HAVE_DAOS)
  target_sources(sst PRIVATE dp/daos_dp.c)
  target_link_libraries(sst PRIVATE DAOS::DAOS)
  set(CMAKE_REQUIRED_INCLUDES ${DAOS_INCLUDE_DIRS})
endif()

if(ADIOS2_HAVE_ZFP)
  target_sources(sst PRIVATE cp/ffs_zfp.c)
  target_link_libraries(sst PRIVATE zfp::zfp)
endif()

if(ADIOS2_SST_BUILD_MPI_DP)
  target_sources(sst PRIVATE dp/mpi_dp.c)
  target_link_libraries(sst PRIVATE MPI::MPI_C)
endif()

# Set library version information
set_target_properties(sst PROPERTIES
  OUTPUT_NAME adios2${ADIOS2_LIBRARY_SUFFIX}_sst
  VERSION ${ADIOS2_LIBRARY_VERSION}
  SOVERSION ${ADIOS2_LIBRARY_SOVERSION}
)

#------------------------------------------------------------------------------#
# Generate SST defines
#------------------------------------------------------------------------------#
set(SST_CONFIG_OPTS
  LIBFABRIC
  UCX
  FI_GNI
  CRAY_DRC
  CRAY_CXI
  NVStream
  MPI
)
include(SSTFunctions)
GenerateSSTHeaderConfig(${SST_CONFIG_OPTS})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SSTConfig.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/adios2/toolkit/sst
  COMPONENT adios2_core-development
)

if (NOT MSVC)
  target_compile_features(sst PRIVATE ${ADIOS2_C99_FEATURES})
endif()

target_include_directories(sst PRIVATE
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/source>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
target_link_libraries(sst PRIVATE adios2::thirdparty::EVPath)

if(NOT BUILD_SHARED_LIBS)
  # When sst is included in a static library we must export the
  # target so that consumers know the link dependencies.
  install(TARGETS sst EXPORT adios2Exports)
endif()

if(CMAKE_VERSION VERSION_LESS 3.14)
  # CMake < 3.14 forgets to propagate private dependencies of object libraries.
  # Propagate them the same way CMake >= 3.14 do.
  # FIXME: Drop this workaround when we require CMake 3.14.
  get_property(sst_deps TARGET sst PROPERTY LINK_LIBRARIES)
  foreach(dep ${sst_deps})
    target_link_libraries(sst INTERFACE "$<LINK_ONLY:${dep}>")
  endforeach()
endif()

if (NOT ADIOS2_USE_PIP)
  add_subdirectory(util)
endif ()
