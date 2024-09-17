#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

if(CMAKE_C_COMPILER_ID MATCHES "^(GNU)$")
  # Silence a noisy warning in our log functions as the behavior is
  # intentional
  string(APPEND CMAKE_C_FLAGS " -Wno-format-zero-length")
endif()

include(CheckIncludeFiles)
check_include_files(unistd.h HAVE_UNISTD_H)
check_include_files(windows.h HAVE_WINDOWS_H)

if (HAVE_UNISTD_H)
   add_definitions(-DHAVE_UNISTD_H)
endif()

if (HAVE_WINDOWS_H)
   add_definitions(-DHAVE_WINDOWS_H)
endif()

add_library(adios2_h5vol
  H5VolReadWrite.c
  H5Vol.c
  H5VolUtil.c
  H5Vol_attr.c
  H5Vol_dataset.c
  H5Vol_file.c
  H5Vol_group.c
  H5Vol_link.c
  H5Vol_object.c
)
target_include_directories(adios2_h5vol PRIVATE ${HDF5_C_INCLUDE_DIRS})	 
target_link_libraries(adios2_h5vol PRIVATE
  ${HDF5_C_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT}
)

if(ADIOS2_USE_MPI)
  target_link_libraries(adios2_h5vol PRIVATE adios2::c_mpi MPI::MPI_C)
else()
  target_link_libraries(adios2_h5vol PRIVATE adios2::c)
endif()

# VOL install should only contain the .so and no headers
install(TARGETS adios2_h5vol
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  COMPONENT adios2_h5vol
)
