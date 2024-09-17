#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

if(ADIOS2_HAVE_MPI AND NOT WIN32) # Doesn't reliably work on windows
  gtest_add_tests_helper(ManyVars MPI_ONLY "" Performance. "")

  # Pure C code, not added to test, 
  # just for executing manually for performance studies
  add_executable(PerfManyVars PerfManyVars.c)
  target_link_libraries(PerfManyVars adios2::c_mpi MPI::MPI_C)
endif()
