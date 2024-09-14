#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

if(ADIOS2_HAVE_MPI)
  add_subdirectory(globalArray1D)
endif()
add_subdirectory(globalArrayND)
add_subdirectory(joinedArray)
add_subdirectory(localArray)
add_subdirectory(queryWorker)
add_subdirectory(values)
add_subdirectory(variablesShapes)
