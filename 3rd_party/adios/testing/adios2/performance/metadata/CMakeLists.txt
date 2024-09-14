#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

if(ADIOS2_HAVE_MPI)
  # just for executing manually for performance studies
  add_executable(PerfMetaData PerfMetaData.cpp)
  target_link_libraries(PerfMetaData adios2::cxx11_mpi MPI::MPI_C)
endif()
