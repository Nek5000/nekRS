#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#

add_subdirectory(engine)
if(ADIOS2_HAVE_Sodium)
  add_subdirectory(operator)
endif()
