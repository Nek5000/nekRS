#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

if(ADIOS2_HAVE_Python)
  add_subdirectory(Python)
endif()

if(ADIOS2_HAVE_Fortran)
  add_subdirectory(Fortran)
endif()

add_subdirectory(C)
add_subdirectory(CXX11)
