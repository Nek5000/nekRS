#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/iotest-config
  DESTINATION ${PROJECT_BINARY_DIR}
)

add_executable(adios_iotest settings.cpp decomp.cpp processConfig.cpp ioGroup.cpp stream.cpp adiosStream.cpp adios_iotest.cpp)
target_link_libraries(adios_iotest adios2::cxx11_mpi MPI::MPI_C adios2_core_mpi)
if(WIN32)
  target_link_libraries(adios_iotest getopt)
endif()
set_property(TARGET adios_iotest PROPERTY OUTPUT_NAME adios2_iotest${ADIOS2_EXECUTABLE_SUFFIX})

if(ADIOS2_HAVE_HDF5 AND HDF5_IS_PARALLEL)
  target_compile_definitions(adios_iotest PRIVATE ADIOS2_HAVE_HDF5_PARALLEL)

  if(HDF5_C_INCLUDE_DIRS)
    target_include_directories(adios_iotest PRIVATE ${HDF5_C_INCLUDE_DIRS})
  else()
    target_include_directories(adios_iotest PRIVATE ${HDF5_INCLUDE_DIRS})
  endif()

  target_sources(adios_iotest PRIVATE hdf5Stream.cpp)

  target_link_libraries(adios_iotest ${HDF5_C_LIBRARIES})

endif()


install(TARGETS adios_iotest EXPORT adios2
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT adios2_tools-runtime
)

install(DIRECTORY iotest-config/
   DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/iotest-config COMPONENT adios2_tools-data
)
