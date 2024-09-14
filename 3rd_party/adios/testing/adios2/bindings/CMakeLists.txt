#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

if(ADIOS2_HAVE_Python)
  add_subdirectory(python)
endif()

if(ADIOS2_HAVE_Fortran)
  add_subdirectory(fortran)
endif()

add_subdirectory(C)
